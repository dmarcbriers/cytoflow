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

'''
Created on Apr 18, 2015

@author: brian
'''

from Queue import PriorityQueue
import heapq

class UniquePriorityQueue(PriorityQueue):
    """
    A PriorityQueue that only allows one copy of each item.
    http://stackoverflow.com/questions/5997189/how-can-i-make-a-unique-value-priority-queue-in-python
    """
    
    def _init(self, maxsize):
        PriorityQueue._init(self, maxsize)
        self.values = set()

    def _put(self, item, heappush=heapq.heappush):
        if item[1] not in self.values:
            self.values.add(item[1])
            PriorityQueue._put(self, item, heappush)
        else:
            pass

    def _get(self, heappop=heapq.heappop):
        item = PriorityQueue._get(self, heappop)
        self.values.remove(item[1])
        return item