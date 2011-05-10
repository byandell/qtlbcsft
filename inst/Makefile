bcsft.html: bcsft.asciidoc
	echo 'library(ascii); Asciidoc("bcsft.asciidoc", width=9, height=4, prefix.string="Figs/fig", driver=cacheSweaveAsciidoc(), pdf=FALSE)' | R --no-save
	asciidoc -a data-uri -a toc -a toclevels=4 -a num bcsft.asciidoc.txt
	mv -f bcsft.asciidoc.html bcsft.html

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt