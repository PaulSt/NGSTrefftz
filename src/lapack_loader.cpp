#include "lapack_loader.hpp"

#include <cstring>
#include <cwchar>

#if defined(_WIN32)
#include <windows.h>
#include <tlhelp32.h>
#elif defined(__APPLE__)
#include <dlfcn.h>
#include <mach-o/dyld.h>
#else
#include <dlfcn.h>
#include <link.h>
#endif

namespace
{
  using RawSymbol = void (*) ();

  RawSymbol SelectPreferredSymbol (RawSymbol ngbla_symbol,
                                   RawSymbol lapack_backend_symbol,
                                   RawSymbol default_symbol,
                                   RawSymbol generic_symbol) noexcept
  {
    if (ngbla_symbol)
      return ngbla_symbol;
    else if (lapack_backend_symbol)
      return lapack_backend_symbol;
    else if (default_symbol)
      return default_symbol;
    else if (generic_symbol)
      return generic_symbol;
    else
      return nullptr;
  }

#if defined(_WIN32)
  using ImageHandle = HMODULE;

  RawSymbol FindSymbolInImage (ImageHandle image, const char *const *names,
                               size_t num_names) noexcept
  {
    for (size_t index = 0; index < num_names; index++)
      {
        FARPROC symbol = GetProcAddress (image, names[index]);
        if (symbol)
          return reinterpret_cast<RawSymbol> (symbol);
      }
    return nullptr;
  }
#else
  bool IsLapackBackend (const char *image) noexcept
  {
    return strstr (image, "mkl") != nullptr
           || strstr (image, "openblas") != nullptr
           || strstr (image, "lapack") != nullptr
           || strstr (image, "accelerate") != nullptr
           || strstr (image, "Accelerate") != nullptr;
  }

  using ImageHandle = void *;

  RawSymbol FindSymbolInImage (ImageHandle image, const char *const *names,
                               size_t num_names) noexcept
  {
    for (size_t index = 0; index < num_names; index++)
      {
        RawSymbol symbol
            = reinterpret_cast<RawSymbol> (dlsym (image, names[index]));
        if (symbol)
          return symbol;
      }
    return nullptr;
  }
#endif

#if defined(_WIN32)
  RawSymbol
  FindLoadedSymbol (const char *const *names, size_t num_names) noexcept
  {
    HANDLE snapshot = CreateToolhelp32Snapshot (
        TH32CS_SNAPMODULE | TH32CS_SNAPMODULE32, GetCurrentProcessId ());
    if (snapshot == INVALID_HANDLE_VALUE)
      return nullptr;

    MODULEENTRY32W module;
    module.dwSize = sizeof (module);
    RawSymbol ngbla_symbol = nullptr;
    RawSymbol lapack_backend_symbol = nullptr;
    RawSymbol generic_symbol = nullptr;
    if (Module32FirstW (snapshot, &module))
      do
        {
          const bool is_ngbla = wcsstr (module.szModule, L"ngbla") != nullptr;
          const bool is_lapack_backend
              = wcsstr (module.szModule, L"mkl") != nullptr
                || wcsstr (module.szModule, L"openblas") != nullptr
                || wcsstr (module.szModule, L"lapack") != nullptr
                || wcsstr (module.szModule, L"accelerate") != nullptr;
          RawSymbol symbol
              = FindSymbolInImage (module.hModule, names, num_names);
          if (is_ngbla && symbol)
            {
              ngbla_symbol = symbol;
              break;
            }
          else if (is_lapack_backend && !lapack_backend_symbol && symbol)
            {
              lapack_backend_symbol = symbol;
            }
          else if (!generic_symbol && symbol)
            {
              generic_symbol = symbol;
            }
        }
      while (Module32NextW (snapshot, &module));
    CloseHandle (snapshot);
    return SelectPreferredSymbol (ngbla_symbol, lapack_backend_symbol, nullptr,
                                  generic_symbol);
  }
#elif defined(__APPLE__)
  RawSymbol
  FindLoadedSymbol (const char *const *names, size_t num_names) noexcept
  {
    RawSymbol ngbla_symbol = nullptr;
    RawSymbol lapack_backend_symbol = nullptr;
    RawSymbol generic_symbol = nullptr;
    for (uint32_t index = 0; index < _dyld_image_count (); index++)
      {
        const char *image = _dyld_get_image_name (index);
        if (!image)
          continue;
        void *handle = dlopen (image, RTLD_LAZY | RTLD_NOLOAD);
        if (!handle)
          continue;
        RawSymbol symbol = FindSymbolInImage (handle, names, num_names);
        dlclose (handle);
        if (strstr (image, "ngbla") && symbol)
          {
            ngbla_symbol = symbol;
            break;
          }
        else if (IsLapackBackend (image) && !lapack_backend_symbol && symbol)
          lapack_backend_symbol = symbol;
        else if (!generic_symbol && symbol)
          generic_symbol = symbol;
      }
    RawSymbol default_symbol = nullptr;
    if (!ngbla_symbol && !lapack_backend_symbol)
      default_symbol = FindSymbolInImage (RTLD_DEFAULT, names, num_names);
    return SelectPreferredSymbol (ngbla_symbol, lapack_backend_symbol,
                                  default_symbol, generic_symbol);
  }
#else
  RawSymbol
  FindLoadedSymbol (const char *const *names, size_t num_names) noexcept
  {
    struct SymbolLookup
    {
      const char *const *names;
      size_t num_names;
      RawSymbol ngbla_symbol = nullptr;
      RawSymbol lapack_backend_symbol = nullptr;
      RawSymbol generic_symbol = nullptr;
    } lookup{ names, num_names };

    dl_iterate_phdr (
        [] (dl_phdr_info *info, size_t, void *data) noexcept -> int {
          // `data` is this C callback's manual capture, passed as &lookup.
          SymbolLookup &lookup = *static_cast<SymbolLookup *> (data);
          if (!info->dlpi_name || !info->dlpi_name[0])
            return 0;

          void *handle = dlopen (info->dlpi_name, RTLD_LAZY | RTLD_NOLOAD);
          if (!handle)
            return 0;
          RawSymbol symbol
              = FindSymbolInImage (handle, lookup.names, lookup.num_names);
          dlclose (handle);
          if (strstr (info->dlpi_name, "ngbla") && symbol)
            {
              lookup.ngbla_symbol = symbol;
              return 1;
            }
          else if (IsLapackBackend (info->dlpi_name)
                   && !lookup.lapack_backend_symbol && symbol)
            lookup.lapack_backend_symbol = symbol;
          else if (!lookup.generic_symbol && symbol)
            lookup.generic_symbol = symbol;
          return 0;
        },
        &lookup);

    RawSymbol default_symbol = nullptr;
    if (!lookup.ngbla_symbol && !lookup.lapack_backend_symbol)
      default_symbol = FindSymbolInImage (RTLD_DEFAULT, names, num_names);
    return SelectPreferredSymbol (lookup.ngbla_symbol,
                                  lookup.lapack_backend_symbol, default_symbol,
                                  lookup.generic_symbol);
  }
#endif

  template <typename FUNCTION, size_t N>
  FUNCTION FindLapackSymbol (const char *const (&names)[N]) noexcept
  {
    RawSymbol symbol = FindLoadedSymbol (names, N);
    if (symbol)
      return reinterpret_cast<FUNCTION> (symbol);
    return nullptr;
  }
}

namespace ngstrefftz
{
  const LapackGesdd &GetLapackGesdd ()
  {
    static constexpr const char *real_names[]
        = { "dgesdd_", "dgesdd", "DGESDD_", "DGESDD" };
    static constexpr const char *complex_names[]
        = { "zgesdd_", "zgesdd", "ZGESDD_", "ZGESDD" };
    static const LapackGesdd symbols{ FindLapackSymbol<RealGesdd> (real_names),
                                      FindLapackSymbol<ComplexGesdd> (
                                          complex_names) };
    return symbols;
  }
}
