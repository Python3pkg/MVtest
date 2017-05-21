import time
import sys

__copyright__ = "Copyright (C) 2015 Todd Edwards, Chun Li and Eric Torstenson"
__license__ = "GPL3.0"
#     This file is part of MVtest.
#
#     MVtest is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     MVtest is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with MVtest.  If not, see <http://www.gnu.org/licenses/>.


class SimpleTimer:
    """Simple abstraction to allow for basic timing. """

    def __init__(self):
        self.start = time.time()

    def report(self, msg, do_reset=False, file=sys.stdout):
        """Print to stdout msg followed by the runtime.

        When true, do_reset will result in a reset of start time.

        """
        print("%s (%s s)" % (msg, time.time() - self.start), file=file)
        if do_reset:
            self.start = time.time()

    def result(self, msg, do_reset=False):
        """Return log message containing ellapsed time as a string.

        When true, do_reset will result in a reset of start time.

        """
        result = "%s (%s s)" % (msg, time.time() - self.start)
        if do_reset:
            self.start = time.time()
        return result

    def reset(self):
        """Reset start time"""
        self.start = time.time()


    def runtime(self):
        """Return ellapsed time and reset start. """
        t = time.time() - self.start
        self.start = time.time()
        return t
