# Ocean Health Index for Canada [can]

Build status of branches:

- [**draft**](https://github.com/OHI-Science/can/tree/draft)

  [![](https://api.travis-ci.org/OHI-Science/can.svg?branch=draft)](https://travis-ci.org/OHI-Science/can/branches)

- [**published**](https://github.com/OHI-Science/can/tree/published)

  not applicable (since merely a copy of a passing draft branch)  

- [**gh-pages**](https://github.com/OHI-Science/can/tree/gh-pages)

  [![](https://api.travis-ci.org/OHI-Science/can.svg?branch=gh-pages)](https://travis-ci.org/OHI-Science/can/branches)
  
  Note that a "build failing" for this gh-pages branch (ie the website) could simply mean a broken link was found.

- [**gh-pages**](https://github.com/OHI-Science/can/tree/app)

  not applicable (because deployment of the Shiny app is done by OHI-Science internally)

For more details, see below and http://ohi-science.org/can/docs.

## gh-pages: website

A "build failing" for the gh-pages branch could simply mean broken links were found.

To test the website locally, install [jekyll](http://jekyllrb.com/docs/installation/) and run:

```bash
jekyll serve --baseurl ''
```

To test links, install html-proofer (`sudo gem install html-proofer`) and run:

```bash
jekyll serve --baseurl ''
# Ctrl-C to stop
htmlproof ./_site
```
