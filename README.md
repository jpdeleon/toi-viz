# toi-viz
TOI visualization plots

## Source
* data: [raw](https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv), [table](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=TOI)
* [code](github.com/jpdeleon/toi-viz)

## Output
* https://jpdeleon.github.io/toi-viz

## Rationale
- TOI list is growing quickly and updated frequently
- There is a need to easily identify the most interesting targets, thru
  - interactive pre-defined plots
  - simple filtering
  - etc

## Framework options
- altair, panel, hvplot
- pyscript
  - pro: no need for webserver; (python) code integrated in html code 
  - con: html file does not render properly in local browser, except in gh-pages
- bokeh
- flask
- ~~streamlit~~ does not output static html

