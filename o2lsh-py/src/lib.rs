#![feature(proc_macro, specialization, const_fn)]
#![feature(const_fn, const_align_of, const_size_of, const_ptr_null, const_ptr_null_mut)]
extern crate pyo3;
use pyo3::prelude::*;

#[py::class]
struct MyClass {
}

#[py::methods]
impl MyClass {

   #[new]
   fn __new__(obj: &PyRawObject) -> PyResult<()> {
       obj.init(|_x| MyClass {})
   }

   fn hello(&self, name: String) -> PyResult<String> {
       Ok(name)
   }
}

#[py::modinit(_o2lsh)]
fn init_mod(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<MyClass>()?;

    Ok(())
}
