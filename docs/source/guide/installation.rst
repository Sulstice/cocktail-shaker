.. _install:

Installation
============

cocktail_shaker supports Python versions 3.3+. There required dependencies, most notably RDKit must be installed.

Option 1: Use pip (recommended)
-------------------------------

The easiest and recommended way to install is using pip::

    pip install cocktail_shaker

This will download the latest version of cocktail_shaker, and place it in your `site-packages` folder so it is automatically
available to all your python scripts.

If you don't already have pip installed, you can `install it using get-pip.py`_::

       curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
       python get-pip.py

Option 2: Download the latest release
-------------------------------------

Alternatively, `download the latest release`_ manually and install yourself::

    tar -xzvf cocktail_shaker-1.0.0tar.gz
    cd cocktail_shaker-1.0.0
    python setup.py install

The setup.py command will install cocktail_shaker in your `site-packages` folder so it is automatically available to all your
python scripts.

Option 3: Clone the repository
------------------------------

The latest development version of cocktail_shaker is always `available on GitHub`_. This version is not guaranteed to be
stable, but may include new features that have not yet been released. Simply clone the repository and install as usual::

    git clone https://github.com/Sulstice/Cocktail-Shaker.git
    cd Cocktail-Shaker
    python setup.py install

.. _`install it using get-pip.py`: http://www.pip-installer.org/en/latest/installing.html
.. _`download the latest release`: https://github.com/mcs07/cocktail_shaker/releases/
.. _`available on GitHub`: https://github.com/Sulstice/Cocktail-Shaker
