#! /bin/bash
# Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#

# This file should contain some kind of test for the octopus
echo ""
echo "=========================="
echo "Running octopus test #2..."

cat<<EOF>inp
%Species
"ho" | 1.0 | usdef | 1 | "0.5*r^2"
%
%Coordinates
"ho" | 0 | 0 | 0 | no
%
spacing = 0.23
radius = 7.0
BoxShape = sphere
PoissonSolver = fft
EOF
../src/oct-test < inp
output=$?

echo "Finished octopus test #2."
echo "========================="
echo ""

exit ${output}