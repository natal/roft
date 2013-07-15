all:
	rust build roft.rc -L./kiss3d/glcore-rs/lib/ -L./kiss3d/glfw-rs/lib/ -L./kiss3d/lib -L./nalgebra/lib -L./nphysics/lib -L./nphysics/ncollide/lib -L./kiss3d/rust-stb-image/ --opt-level=3
	rm -f blob.dot

deps:
	make -C nalgebra
	make deps -C nphysics
	make deps -C kiss3d
	make -C nphysics
	make -C kiss3d

simple:
	dot -Kfdp -n -Tpdf simple.dot -o simple.pdf
	open simple.pdf

line:
	dot -Kfdp -n -Tpdf line.dot -o line.pdf
	open line.pdf

blob:
	dot -Kfdp -n -Tpdf blob.dot -o blob.pdf
	open blob.pdf
