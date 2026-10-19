.PHONY: setup example run test lint clean

VENV ?= ../venv
PY   := $(VENV)/bin/python

setup:
	python3 -m venv $(VENV)
	$(VENV)/bin/pip install -q --upgrade pip
	$(VENV)/bin/pip install -q -r requirements.txt

example:
	$(PY) data/make_example_data.py --outdir data/example

run:
	$(PY) src/pipeline.py --counts data/example/counts.csv \
		--metadata data/example/metadata.csv --outdir results \
		--truth data/example/truth.csv

test:
	$(PY) -m pytest tests -q

clean:
	rm -rf results data/example __pycache__ src/__pycache__ tests/__pycache__
