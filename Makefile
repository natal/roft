libs=-L./lib/kiss3d/glcore-rs/lib/ -L./lib/kiss3d/glfw-rs/lib/ -L./lib/kiss3d/lib -L./lib/nalgebra/lib -L./lib/nphysics/lib -L./lib/nphysics/ncollide/lib -L./lib/kiss3d/rust-stb-image/
libs_w_cl=$(libs) -L./lib/rust-opencl/ -L./lib/rs2cl/lib

all:
	mkdir -p bin
	rust build src/roft_gpu.rc --opt-level=3 $(libs_w_cl) --out-dir bin
	rust build src/roft.rc --opt-level=3 $(libs) --out-dir bin

deps:
	make -C lib/rust-opencl
	make -C lib/nalgebra
	make deps -C lib/nphysics
	make deps -C lib/kiss3d
	make deps -C lib/rs2cl
	make -C lib/nphysics
	make -C lib/kiss3d
	make -C lib/rs2cl
