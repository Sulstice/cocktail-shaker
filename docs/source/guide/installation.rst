.. _install:

Installation
============

cocktail-shaker supports Python versions 3.3+. Pyhton's required dependencies, most notably RDKit, must be installed.

Option 1 (Recommended): Use Pip 
-------------------------------

The easiest way to install cocktail-shaker is using pip::

    pip install cocktail_shaker

This will download the latest version of cocktail-shaker and place it in your `site-packages` folder so it is automatically
available to all your Python scripts.

If you do not have pip installed yet, you can `install it using get-pip.py`_::

       curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
       python get-pip.py

Option 2: Download the Latest Release
-------------------------------------

Alternatively, you can get cocktail-shaker by manually `download the latest release`_ and installing it yourself::

    tar -xzvf cocktail_shaker-1.0.0tar.gz
    cd cocktail_shaker-1.0.0
    python setup.py install

The setup.py command will install cocktail-shaker in your `site-packages` folder so it is automatically available to all your
Python scripts.

Option 3: Clone the Repository
------------------------------

Lastly, the latest development version of cocktail-shaker is always `available on GitHub`_. The version on GitHub is not guaranteed to be stable, but may include new features that have not yet been released. 

Simply clone the repository and install as usual::

    git clone https://github.com/Sulstice/Cocktail-Shaker.git
    cd Cocktail-Shaker
    python setup.py install

.. _`install it using get-pip.py`: http://www.pip-installer.org/en/latest/installing.html
.. _`download the latest release`: https://github.com/mcs07/cocktail_shaker/releases/
.. _`available on GitHub`: https://github.com/Sulstice/Cocktail-Shaker
