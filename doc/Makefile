POD = abcsim.pod files.pod mau2cfg.pod tabcfg.pod \
      abcml.pod  cplcfg.pod intro.pod 


all : txt html latex

latex : dummy
	echo Making latex
	cd pod; for i in $(POD); do pod2latex $$i; done
	mv pod/*.tex latex

txt : dummy
	echo Making txt
	cd pod; for i in $(POD); do echo $$i; \
     perldoc -t $$i > ../txt/`basename $$i .pod`.txt; done

html : dummy
	echo Making html
	cd pod; for i in $(POD); do echo $$i; \
  pod2html --infile=$$i > ../html/`basename $$i .pod`.html; done

clean : dummy
	rm -rf txt/* html/*

dummy :
