.DELETE_ON_ERROR:
GENICE=genice2
BASE=genice2_cage
PACKAGE=genice2-cage

all: README.md

%: temp_% replacer.py $(BASE)/formats/cage.py $(BASE)/__init__.py
	python replacer.py < $< > $@

test: 1h.cage.gro.test 1h.cage.yap.test CS1.cage.test CS2.cage.json.test CRN1.cage.json.test CRN2.cage.test
CS1.cage: $(BASE)/formats/cage.py Makefile
	( cd $(BASE) && $(GENICE) CS1 -r 2 2 2 -f cage[12,14-16:maxring=6] ) > $@
CS2.cage.json: $(BASE)/formats/cage.py Makefile
	( cd $(BASE) && $(GENICE) CS2 -f cage[12,14-16:maxring=6:json] ) > $@
CRN1.cage.json: $(BASE)/formats/cage.py Makefile
	( cd $(BASE) && $(GENICE) CRN1 -f cage[-12:maxring=8:json] ) > $@
CRN2.cage: $(BASE)/formats/cage.py Makefile
	( cd $(BASE) && $(GENICE) CRN2 -f cage[ring=6] ) > $@
1h.cage.yap: $(BASE)/formats/cage.py Makefile
	( cd $(BASE) && $(GENICE) 1h -r 2 2 2 -f cage[-5:maxring=6:yaplot] ) > $@
1h.cage.gro: $(BASE)/formats/cage.py Makefile
	( cd $(BASE) && $(GENICE) 1h -r 2 2 2 -f cage[-5:maxring=6:gromacs] ) > $@
%.test:
	make $*
	diff $* ref/$*


test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install --index-url https://test.pypi.org/simple/ $(PACKAGE)



install:
	./setup.py install
uninstall:
	-pip uninstall -y $(PACKAGE)
build: README.md $(wildcard $(BASE)/formats*.py)
	./setup.py sdist bdist_wheel


deploy: build
	twine upload --repository pypi dist/*
check:
	./setup.py check
clean:
	-rm $(ALL) *~ */*~ *.cage *.json
	-rm -rf build dist *.egg-info
	-find . -name __pycache__ | xargs rm -rf
