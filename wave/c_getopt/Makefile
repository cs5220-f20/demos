wave_demo.x: wave_demo.c
	gcc -O3 -ffast-math wave_demo.c -o wave_demo.x

%.md: %.c
	ldoc wave_demo.c -o wave_demo.md

%.pdf: %.md
	pandoc $< -o $@

%.html: %.md
	pandoc --mathjax --standalone $< -o $@

clean:
	rm -f wave_demo.x
