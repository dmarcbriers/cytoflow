#!/usr/bin/env python2.7

# (c) Massachusetts Institute of Technology 2015-2016
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import absolute_import

import matplotlib as mpl
import seaborn as sns

mpl.rc('legend', markerscale = 5)
sns.set_style("whitegrid", {
                "xtick.major.size": 6,
                "ytick.major.size": 6,
                "xtick.minor.size": 3,
                "ytick.minor.size": 3,
                })
sns.set_context("talk")

from .i_view import IView
from .i_selectionview import ISelectionView

from .bar_chart import BarChartView
from .hexbin import HexbinView
from .histogram import HistogramView
from .scatterplot import ScatterplotView
from .stats_1d import Stats1DView
from .stats_2d import Stats2DView