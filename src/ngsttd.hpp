#ifndef NGSTTD
#define NGSTTD

#ifdef WIN32
#define NGSTREFFTZ_API_EXPORT __declspec (dllexport)
#define NGSTREFFTZ_API_IMPORT __declspec (dllimport)
#else
#define NGSTREFFTZ_API_EXPORT
#define NGSTREFFTZ_API_IMPORT
#endif

#ifdef NGSTREFFTZ_EXPORTS
#define NGST_DLL NGSTREFFTZ_API_EXPORT
#else
#define NGST_DLL NGSTREFFTZ_API_IMPORT
#endif

#endif
