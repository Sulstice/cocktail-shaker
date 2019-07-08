
## Ligand Library Design
![Build status](https://travis-ci.com/Sulstice/Ligand-Library-Creation.svg?branch=master)
![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)


**Background:** ML algorithms are capable of searching extensive chemical spaces
for the best starting points for a drug development program. This leaves us with the question of
how do we design and construct the virtual library which will be searched during this process.

**Goal:** Construct a library of 10,000 distinct compounds, which are all synthesizable. The library
should be optimized for chemical diversity, to span as much of the chemical space as possible.

This problem target the main apex pain point hindering chemical modeling today. Given the rise of the genomic and protenoimics have discovered therapeutic targets with no small molecule modulators. The demand asks for more increased virtual high throughput screening but with a lack of adequate screening libraries. 

Questions arise:
    - How do we construct these chemical libraries?
    
    - How do we construct these chemical libraries expanding the most chemical diversity possible?
    
    - How do we construct chemical libraries exploring as many combinations of representation of data as possible? 
    
    - Added on to this what style of chemical libaries are we constructing? (Natural/Non-Natural Amino Acids versus small molecule inhibitors) 
    
    - How do we create these chemical libraries so it adheres to easy synthesis or, to make the chemists happy, can be constructed through click chemistry methods.
    
    - How do construct datasets that are not in violation with existing patents for pharma companies today? 

#### Ligand Selection


