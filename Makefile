PKGNAME=genice_cage
.DELETE_ON_ERROR:

all: README.md

%: temp_% replacer.py $(PKGNAME)/formats/cage.py $(PKGNAME)/__init__.py
	python replacer.py < $< > $@

test: 1h.cage.yap.test CS1.cage.test CS2.cage.json.test CRN1.cage.json.test 
CS1.cage: $(PKGNAME)/formats/cage.py Makefile
	genice CS1 -r 2 2 2 -f cage[12,14-16:maxring=6] > $@
CS2.cage.json: $(PKGNAME)/formats/cage.py Makefile
	genice CS2 -f cage[12,14-16:maxring=6:json] > $@
CRN1.cage.json: $(PKGNAME)/formats/cage.py Makefile
	genice CRN1 -f cage[-12:maxring=8:json] > $@
1h.cage.yap: $(PKGNAME)/formats/cage.py Makefile
	genice 1h -r 2 2 2 -f cage[-5:maxring=6:yaplot] > $@
%.test:
	make $*
	diff $* ref/$*


test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install pillow
	pip install --index-url https://test.pypi.org/simple/ $(PKGNAME)



install:
	./setup.py install
uninstall:
	-pip uninstall -y genice-cage
build: README.md $(wildcard $(PKGNAME)/formats*.py)
	./setup.py sdist bdist_wheel


deploy: build
	twine upload --repository pypi dist/*
check:
	./setup.py check
clean:
	-rm $(ALL) *~ */*~ *.cage *.json
	-rm -rf build dist *.egg-info
	-find . -name __pycache__ | xargs rm -rf
