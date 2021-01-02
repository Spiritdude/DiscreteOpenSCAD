NAME = DiscreteOpenSCAD
VERSION = 0.0.6

all::
	@echo "make edit install deinstall (init change push pull backup)"

edit::
	dee4 discrete.scad example*.scad single.scad Makefile

install::
	mkdir ~/lib/openscad
	cp -p discrete.scad ~/lib/openscad/

deinstall::
	rm -f ~/lib/openscad/discrete.scad

# -- devs only

init::
	git remote add origin git@github.com:Spiritdude/DiscreteOpenSCAD.git

change::
	git commit -am "..."

push::
	git branch -M main
	git push -u origin main

pull::
	git pull
   
backup::
	cd ..; tar cfz ~/Backup/${NAME}-${VERSION}.tar.gz ${NAME}; scp ~/Backup/${NAME}-${VERSION}.tar.gz backup:Backup/

