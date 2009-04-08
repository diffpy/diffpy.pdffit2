OUTPUT = ../build/epydoc-html
WEBDOCDIR = login.cacr.caltech.edu:/cacr/home/proj/danse/docs.danse.us/docroot/diffraction

all:
	epydoc --config=epydoc.cfg --output=$(OUTPUT)

clean:
	/bin/rm -rf $(OUTPUT)


publish:
	rsync -e /usr/bin/ssh -pRt --chmod=a+rwx,o-w,Dg+s,Fa-x \
	    $(OUTPUT)/ $(WEBDOCDIR)/diffpy.pdffit2
