.DELETE_ON_ERROR:
GENICE=genice2
BASENAME=genice2_cage
PIPNAME=genice2-cage

all: README.md

test: 1h.cage.gro.test 1h.cage.yap.test CS1.cage.test CS2.cage.json.test CRN1.cage.json.test CRN2.cage.test zra-d.solid.json.test
CS1.cage: $(BASENAME)/formats/cage.py Makefile
	( cd $(BASENAME) && $(GENICE) CS1 -r 2 2 2 -f cage[12,14-16:maxring=6] ) > $@
CS2.cage.json: $(BASENAME)/formats/cage.py Makefile
	( cd $(BASENAME) && $(GENICE) CS2 -f cage[12,14-16:maxring=6:geom] ) > $@
CRN1.cage.json: $(BASENAME)/formats/cage.py Makefile
	( cd $(BASENAME) && $(GENICE) CRN1 -f cage[-12:maxring=8:geom] ) > $@
CRN2.cage: $(BASENAME)/formats/cage.py Makefile
	( cd $(BASENAME) && $(GENICE) CRN2 -f cage[ring=6] ) > $@
1h.cage.yap: $(BASENAME)/formats/cage.py Makefile
	( cd $(BASENAME) && $(GENICE) ice1h -r 2 2 2 -f cage[-5:maxring=6:yaplot] ) > $@
1h.cage.gro: $(BASENAME)/formats/cage.py Makefile
	( cd $(BASENAME) && $(GENICE) ice1h -r 2 2 2 -f cage[-5:maxring=6:gromacs] ) > $@
zra-d.solid.json: $(BASENAME)/formats/cage.py Makefile
	( cd $(BASENAME) && $(GENICE) zra-d -r 2 2 2 -f cage[solid] ) > $@
FAU.cage.yap: $(BASENAME)/formats/cage.py Makefile
	( cd $(BASENAME) && $(GENICE) FAU -r 2 2 2 -f cage[-30:maxring=6:yaplot] ) > $@
%.test:
	make $*
	diff $* ref/$*


%: temp_% replacer.py $(wildcard $(BASENAME)/formats/*.py) $(wildcard $(BASENAME)/*.py) pyproject.toml
	python replacer.py $< > $@


test-deploy: clean
	poetry publish --build -r testpypi
test-install:
	pip install --index-url https://test.pypi.org/simple/ $(PIPNAME)
uninstall:
	-pip uninstall -y $(PIPNAME)
build: README.md $(wildcard genice2_yaplot/*.py)
	poetry build
deploy: clean
	poetry publish --build
check:
	poetry check


clean:
	-rm $(ALL) *~ */*~ *svg CS2.png
	-rm -rf build dist *.egg-info
	-find . -name __pycache__ | xargs rm -rf

distclean: clean
	rm -rf *.json *.gro *.yap *.cage @*
