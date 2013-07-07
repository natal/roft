all:
	rustc -L./nalgebra/lib graph.rs --opt-level=3
	rm -rf out.dot
	./graph
	dot -Kfdp -n -Tpdf out.dot -o outfile.pdf
	open outfile.pdf
