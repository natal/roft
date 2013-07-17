libs=-L./kiss3d/glcore-rs/lib/ -L./kiss3d/glfw-rs/lib/ -L./kiss3d/lib -L./nalgebra/lib -L./nphysics/lib -L./nphysics/ncollide/lib -L./kiss3d/rust-stb-image/
libs_w_cl=$(libs) -L./rust-opencl/ -L./rs2cl/lib

all:
	rust build roft_gpu.rc --opt-level=3 $(libs_w_cl)
	rust build roft.rc --opt-level=3 $(libs)

deps:
	make -C rust-opencl
	make -C nalgebra
	make deps -C nphysics
	make deps -C kiss3d
	make -C nphysics
	make -C kiss3d
	make -C rs2cl

simple:
	dot -Kfdp -n -Tpdf simple.dot -o simple.pdf
	open simple.pdf

line:
	dot -Kfdp -n -Tpdf line.dot -o line.pdf
	open line.pdf

blob:
	dot -Kfdp -n -Tpdf blob.dot -o blob.pdf
	open blob.pdf
