#include <iostream>
#include <cmath>
#include <array>
#include <cstdint>

using namespace std;

#include "MultiArray.hpp"


namespace ngfem
{
	template <class T, int depth>
	void MultiArray<T, depth> :: put(array<int, depth> indices, T thing_to_put)
	{
		int index = 0;
		int dimension_counter = 1;
		for (int i : indices)
		{
			if(i >= each_size)
			{
				cout << "error" << endl;
				return;
			}
			else
			{
				index += i * dimension_counter;
				dimension_counter *= each_size;
			}
		}
		this->multiArray[index] = thing_to_put;
		//cout << "\n At index " << index << " I put: " << this->multiArray[index] << endl;
		//cout << "\n At index " << indices[0] << " " << indices[1] << " " << indices[2] << " " << " I put: " << this->multiArray[index] << endl;
	}

	template <class T, int depth>
	T MultiArray<T, depth> :: get(array<int, depth> indices) const
	{
		int index = 0;
		int dimension_counter = 1;
		for (int i : indices)
		{
			if(i >= each_size)
			{
				cout << "error" << endl;
			}
			else
			{
				index += i * dimension_counter;
				dimension_counter *= each_size;
			}
		}
		T thing_to_get = this->multiArray[index];
		//cout << "\n At index "<< index << " I found: " << thing_to_get << endl;
		return thing_to_get;
	}

	template <class T, int depth>
	T MultiArray<T, depth> :: operator[](const std::initializer_list<int> indices)
	{
		int index = 0;
		int dimension_counter = 1;
		for (int i : indices)
		{
			if(i >= each_size)
			{
				cout << "error" << endl;
				break;
			}
			else
			{
				index += i * dimension_counter;
				dimension_counter *= each_size;
			}
		}
		return this->multiArray[index];
	}

}



#ifdef NGS_PYTHON
void ExportMultiArray(py::module m)
{
  using namespace ngfem;
  py::class_<MultiArray<float, 3>, shared_ptr<MultiArray<float, 3> > >
    (m, "MultiArray", "test")
    //.def(py::init<>())
		.def(py::init<int>())
		.def("put", &MultiArray<float, 3>::put)
		.def("get", &MultiArray<float, 3>::get)
		;
}
#endif // NGS_PYTHON
