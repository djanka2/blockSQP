blockSQP - a block structure exploiting SQP solver
Copyright (c) 2012-2015 Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>

Introduction
============


Documentation
=============


Getting blockSQP
================



Building and Installation
=========================

* Get ma57 from HSL

* add `-fPIC` to the FFLAGS and FCFLAGS in Makefile and SRC/Makefile after
  executing configure.

* Install ma57 using make and make install

* Get special Schur complement branch from qpOASES::

  svn co https://projects.coin-or.org/svn/qpOASES/branches/schurComplement/ .

* Compile qpOASES using make

* Add qpOASES directory to .bashrc or to Makefile as variable `QPOASESDIR`.

* Compile and install blocksqp using make


Licensing
=========

The library is published under the very permissive zlib free software
license which should allow you to use the software wherever you need.

This is the full license text (zlib license):

    blockSQP
    Copyright (c) 2012-2015 Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>

    This software is provided 'as-is', without any express or implied
    warranty. In no event will the authors be held liable for any damages
    arising from the use of this software.

    Permission is granted to anyone to use this software for any purpose,
    including commercial applications, and to alter it and redistribute it
    freely, subject to the following restrictions:

       1. The origin of this software must not be misrepresented; you must not
       claim that you wrote the original software. If you use this software
       in a product, an acknowledgment in the product documentation would be
       appreciated but is not required.

       2. Altered source versions must be plainly marked as such, and must not
       be misrepresented as being the original software.

       3. This notice may not be removed or altered from any source
       distribution.


Acknowledgements
================

