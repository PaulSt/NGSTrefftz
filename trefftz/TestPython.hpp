#ifndef FILE_TESTPYTHON_HPP
#define FILE_TESTPYTHON_HPP
#include <Python.h>
#include <python_ngstd.hpp>
#include <pybind11/iostream.h>

int TestPython ()
{
  PyObject *pName, *pModule, *pDict, *pFunc;
  py::module test = py::module::import ("EvolveTent");
  py::object bla = test.attr ("add") (1, 2);
  // Py_Initialize();
  // PyRun_SimpleString("import sys\nsys.path.append(\".\")");
  // const char* plugInCode = "EvolveTent";
  // pName = PyUnicode_FromString(plugInCode);
  // pModule = PyImport_Import(pName);
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
