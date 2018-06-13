#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
#include <Python.h>
#include <python_ngstd.hpp>
#include <pybind11/iostream.h>

int TestPython ()
{
  PyObject *pName, *pModule, *pDict, *pFunc;
  string PlugInPath = "~/trefftzngs/trefftz/EvolveTent.py";
  string plugInInitCode;
  // Py_SetProgramName(argv[0]);  [> optional but recommended <]
  Py_Initialize ();
  PyRun_SimpleString ("from time import time,ctime\n"
                      "print('Today is',ctime(time()))\n");
  plugInInitCode = "import sys\nimport os\nsys.path.insert(0,\"";
  plugInInitCode += PlugInPath;
  plugInInitCode += "\")";
  PyRun_SimpleString ("import sys\nsys.path.append(\".\")");
  const char *plugInCode = "EvolveTent";
  pName = PyUnicode_FromString (plugInCode);
  pModule = PyImport_Import (pName);
  return 4;
}

#ifdef NGS_PYTHON
void ExportTestPython (py::module m)
{
  m.def ("TestPythonF",
         [] () -> int {
           return TestPython ();
         } //, py::call_guard<py::gil_scoped_release>()
  );
};
#endif // NGS_PYTHON

#endif
