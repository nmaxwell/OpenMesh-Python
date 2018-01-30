#ifndef OPENMESH_PYTHON_VECTOR_HH
#define OPENMESH_PYTHON_VECTOR_HH

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


template <class Vector, class Scalar>
void set_item(Vector& _vec, int _index, Scalar _value) {
	if (_index < 0) {
		_index += _vec.size();
	}

	if ((size_t)_index < _vec.size()) {
		_vec[_index] = _value;
	}
	else {
		throw py::index_error();
	}
}

template <class Vector, class Scalar>
Scalar get_item(Vector& _vec, int _index) {
	if (_index < 0) {
		_index += _vec.size();
	}

	if ((size_t)_index < _vec.size()) {
		return _vec[_index];
	}
	else {
		throw py::index_error();
	}

	return 0.0;
}

namespace {
template<class Scalar>
struct Factory {
	typedef OpenMesh::VectorT<Scalar, 2> Vector2;
	typedef OpenMesh::VectorT<Scalar, 3> Vector3;
	typedef OpenMesh::VectorT<Scalar, 4> Vector4;

	static Vector2 *vec2_default() {
		return new Vector2(Scalar(), Scalar());
	}
	static Vector2 *vec2_user_defined(const Scalar& _v0, const Scalar& _v1) {
		return new Vector2(_v0, _v1);
	}
	static Vector3 *vec3_default() {
		return new Vector3(Scalar(), Scalar(), Scalar());
	}
	static Vector3 *vec3_user_defined(const Scalar& _v0, const Scalar& _v1, const Scalar& _v2) {
		return new Vector3(_v0, _v1, _v2);
	}
	static Vector4 *vec4_default() {
		return new Vector4(Scalar(), Scalar(), Scalar(), Scalar());
	}
	static Vector4 *vec4_user_defined(const Scalar& _v0, const Scalar& _v1, const Scalar& _v2, const Scalar& _v3) {
		return new Vector4(_v0, _v1, _v2, _v3);
	}
};
}

template<class Scalar, class Vector>
void defInitMod(py::module& m, py::class_< OpenMesh::VectorT<Scalar, 2> > &classVector) {
	classVector
		.def(py::init(&Factory<Scalar>::vec2_default))
		.def(py::init(&Factory<Scalar>::vec2_user_defined))
		;
}

template<class Scalar, class Vector>
void defInitMod(py::module& m, py::class_< OpenMesh::VectorT<Scalar, 3> > &classVector) {
	Vector (Vector::*cross)(const Vector&) const = &Vector::operator%;
	classVector
		.def(py::init(&Factory<Scalar>::vec3_default))
		.def(py::init(&Factory<Scalar>::vec3_user_defined))
		.def("__mod__", cross)
		;
	m.def("cross", cross);
}

template<class Scalar, class Vector>
void defInitMod(py::module& m, py::class_< OpenMesh::VectorT<Scalar, 4> > &classVector) {
	classVector
		.def(py::init(&Factory<Scalar>::vec4_default))
		.def(py::init(&Factory<Scalar>::vec4_user_defined))
		;
}

/**
 * Expose a vector type to %Python.
 *
 * This function template is used to expose vectors to %Python. The template
 * parameters are used to instantiate the appropriate vector type.
 *
 * @tparam Scalar A scalar type.
 * @tparam N The dimension of the vector.
 *
 * @param _name The name of the vector type to be exposed.
 *
 * @note N must be either 2, 3 or 4.
 */
template<class Scalar, int N>
void expose_vec(py::module& m, const char *_name) {
	typedef OpenMesh::VectorT<Scalar, N> Vector;

	Scalar (Vector::*min_void)() const = &Vector::min;
	Scalar (Vector::*max_void)() const = &Vector::max;

	Vector (Vector::*max_vector)(const Vector&) const = &Vector::max;
	Vector (Vector::*min_vector)(const Vector&) const = &Vector::min;

	Scalar  (Vector::*dot           )(const Vector&) const = &Vector::operator|;
	Scalar  (Vector::*norm          )(void         ) const = &Vector::norm;
	Scalar  (Vector::*length        )(void         ) const = &Vector::length;
	Scalar  (Vector::*sqrnorm       )(void         ) const = &Vector::sqrnorm;
	Vector& (Vector::*normalize     )(void         )       = &Vector::normalize;
	Vector& (Vector::*normalize_cond)(void         )       = &Vector::normalize_cond;

	Vector& (Vector::*op_selfmul_scalar)(const Scalar&)       = &Vector::operator*=;
	Vector& (Vector::*op_selfmul_vector)(const Vector&)       = &Vector::operator*=;
	Vector& (Vector::*op_selfdiv_scalar)(const Scalar&)       = &Vector::operator/=;
	Vector& (Vector::*op_selfdiv_vector)(const Scalar&)       = &Vector::operator/=;
	Vector& (Vector::*op_selfadd_vector)(const Vector&)       = &Vector::operator+=;
	Vector& (Vector::*op_selfsub_vector)(const Vector&)       = &Vector::operator-=;
	Vector  (Vector::*op_mul_scalar    )(const Scalar&) const = &Vector::operator*;
	Vector  (Vector::*op_mul_vector    )(const Vector&) const = &Vector::operator*;
	Vector  (Vector::*op_div_scalar    )(const Scalar&) const = &Vector::operator/;
	Vector  (Vector::*op_div_vector    )(const Vector&) const = &Vector::operator/;
	Vector  (Vector::*op_add_vector    )(const Vector&) const = &Vector::operator+;
	Vector  (Vector::*op_sub_vector    )(const Vector&) const = &Vector::operator-;
	Vector  (Vector::*op_unary_minus   )(void         ) const = &Vector::operator-;

#if (_MSC_VER >= 1800 || __cplusplus > 199711L || defined(__GXX_EXPERIMENTAL_CXX0X__)) && !defined(OPENMESH_VECTOR_LEGACY)
	Vector (Vector::*normalized)() const = &Vector::normalized;
#else
	const Vector (Vector::*normalized)() const = &Vector::normalized;
#endif

	py::class_<Vector> classVector = py::class_<Vector>(m, _name);

	classVector
		.def("__setitem__", &set_item<Vector, Scalar>)
		.def("__getitem__", &get_item<Vector, Scalar>)
		.def("__eq__", &Vector::operator==)
		.def("__ne__", &Vector::operator!=)
		.def("__lt__", &Vector::operator<)
		.def("__imul__", op_selfmul_scalar)
		.def("__imul__", op_selfmul_vector)
		.def("__itruediv__", op_selfdiv_scalar)
		.def("__itruediv__", op_selfdiv_vector)
		.def("__iadd__", op_selfadd_vector)
		.def("__isub__", op_selfsub_vector)
		.def("__mul__", op_mul_scalar)
		.def("__mul__", op_mul_vector)
		.def("__rmul__", op_mul_scalar)
		.def("__truediv__", op_div_scalar)
		.def("__truediv__", op_div_vector)
		.def("__add__", op_add_vector)
		.def("__sub__", op_sub_vector)
		.def("__neg__", op_unary_minus)
		.def("__or__", dot)

		.def("vectorize", &Vector::vectorize, py::return_value_policy::reference_internal)

		.def("dot", dot)
		.def("norm", norm)
		.def("length", length)
		.def("sqrnorm", sqrnorm)
		.def("normalized", normalized)
		.def("normalize", normalize, py::return_value_policy::reference_internal)
		.def("normalize_cond", normalize_cond, py::return_value_policy::reference_internal)

		.def("l1_norm", &Vector::l1_norm)
		.def("l8_norm", &Vector::l8_norm)

		.def("max", max_void)
		.def("max_abs", &Vector::max_abs)
		.def("min", min_void)
		.def("min_abs", &Vector::min_abs)
		.def("mean", &Vector::mean)
		.def("mean_abs", &Vector::mean_abs)
		.def("minimize", &Vector::minimize, py::return_value_policy::reference_internal)
		.def("minimized", &Vector::minimized)
		.def("maximize", &Vector::maximize, py::return_value_policy::reference_internal)
		.def("maximized", &Vector::maximized)
		.def("min", min_vector)
		.def("max", max_vector)

		.def_static("size", &Vector::size)
		.def_static("vectorized", &Vector::vectorized)

		.def("__init__", [](Vector& _self, py::array_t<Scalar, py::array::c_style | py::array::forcecast> _arr) {
				if (_arr.size() != N) {
					throw std::runtime_error("Incompatible array size!");
				}
				new (&_self) Vector(_arr.data());
			});
		;

	defInitMod<Scalar, Vector>(m, classVector);
}

#endif
