def load_ipython_extension(ipython):
    from diffpy.pdffit2 import PdfFit
    pf = PdfFit()
    print('        Type  help(pdffit)  or  help(topic)  for information.')
    print
    ns = dict(pdffit=PdfFit)
    pf._exportAll(ns)
    ipython.user_ns.update(ns)
    return
