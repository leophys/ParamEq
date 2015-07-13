An ubersimple parametric (graphical) equalizer
=====================================

This is the result of a month of experiments.
This very simple tool answers my questions about
the use of FFT in sound processing.

Use
---

Simply open the file with Wolfram Mathematica and execute it.
You will have at disposal a function that outputs a (poor)
graphical interface to manipulate the input waveform (that has to
be recorded or produced at 44.1kH). The syntax is

ParamEqualizer[waveform]

One can use 3 (slightly modified) band pass and a low shelf and a high
shelf. After the manipulation, clicking on the “Get sound!” button, will
be output the processed sound.


Issues
------

- Definitely not real time
- Slow response


Licence
------


This piece of software is distributed under GPL v3.0 Licence or above.
