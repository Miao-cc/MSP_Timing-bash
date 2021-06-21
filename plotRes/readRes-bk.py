"""Tools for working with pulse time-of-arrival (TOA) data.

In particular, single TOAs are represented by :class:`pint.toa.TOA` objects, and if you
want to manage a collection of these we recommend you use a :class:`pint.toa.TOAs` object
as this makes certain operations much more convenient. You probably want to load one with
:func:`pint.toa.get_TOAs`.
"""
import os
import re
import copy
import hashlib
import numpy.ma
import numpy as np
import astropy.units as u
import astropy.time as time
import astropy.table as table
from pint.pulsar_mjd import Time

# read the residual
def get_Res(timfile):
    """Load and prepare TOAs for PINT use.
    """
    t = Res(timfile)
    return t

def build_table(res, filename=None):
    mjds, mjd_floats, res, errors = zip(
        *[
            (
                r[0],                    # mjd
                r[1],                    # mjd float
                r[2].to_value(u.us),     # residual 
                r[3].to_value(u.us),     # error
            )
            for r in res
        ]
    )
    return table.Table(
        [
            np.arange(len(mjds)),
            table.Column(mjds),
            np.array(mjd_floats) * u.d,
            np.array(res) * u.us,
            np.array(errors) * u.us,
        ],
        names=(
            "index",
            "mjd",
            "mjd_float",
            "res",
            "error",
        ),
        meta={"filename": filename},
    ).group_by("obs")



def read_res_file(filename):
    """Read Residual from the res.dat.
    """
    if isinstance(filename, str):
        with open(filename, "r") as f:
            return read_res_file(f)
    else:
        f = filename

    res = []
    for line in f.readlines():
        fields = line.split() 
        mjd, residual, error = fields[0], fields[1], fields[3]
        residual = float(residual) * u.microsecond
        error = float(error) * u.microsecond
        ii, ff = mjd.split('.')
        t = time.Time(ii, jj, scale='utc', format="mjd", precision=9)
        res.append(mjd, t, residual, error)

def _group_by_gaps(t, gap):
    ix = np.argsort(t)
    t_sorted = t[ix]
    gaps = np.diff(t_sorted)
    gap_starts = np.where(gaps >= gap)[0]
    gsi = np.concatenate(([0], gap_starts + 1, [len(t)]))
    groups_sorted = np.repeat(np.arange(len(gap_starts) + 1), np.diff(gsi))
    groups = np.zeros(len(t), dtype=int)
    groups[ix] = groups_sorted
    return groups

class Res:

    def __init__(self, toafile=None, toalist=None):
        # First, just make an empty container
        self.filename = None

        reslist = read_toa_file(resfile)

        # Check to see if there were any INCLUDEs:
        inc_fns = [x[0][1] for x in self.commands if x[0][0].upper() == "INCLUDE"]
        self.filename = [toafile] + inc_fns if inc_fns else toafile

        if toalist is None:
            raise ValueError("No TOAs found!")
        else:
            if not isinstance(toalist, (list, tuple)):
                raise ValueError("Trying to initialize TOAs from a non-list class")

        self.table = build_table(toalist, filename=self.filename)
        groups = self.get_groups()
        self.table.add_column(groups, name="groups")
        # Add pulse number column (if needed) or make PHASE adjustments

        # We don't need this now that we have a table

    def __len__(self):
        return len(self.table)

    def __getitem__(self, index):
        if not hasattr(self, "table"):
            raise ValueError("This TOAs object is incomplete and does not have a table")
        if isinstance(index, np.ndarray) and index.dtype == np.bool:
            r = copy.deepcopy(self)
            r.table = r.table[index]
            if len(r.table) > 0:
                r.table = r.table.group_by("obs")
            return r
        elif (
            isinstance(index, np.ndarray)
            and index.dtype == np.int
            or isinstance(index, list)
        ):
            r = copy.deepcopy(self)
            r.table = r.table[index]
            if len(r.table) > 0:
                r.table = r.table.group_by("obs")
            return r
        elif isinstance(index, slice):
            r = copy.deepcopy(self)
            r.table = r.table[index]
            if len(r.table) > 0:
                r.table = r.table.group_by("obs")
            return r
        elif isinstance(index, int):
            raise ValueError("TOAs do not support extraction of TOA objects (yet?)")
        else:
            raise ValueError("Unable to index TOAs with {}".format(index))

    def __eq__(self, other):
        sd, od = self.__dict__.copy(), other.__dict__.copy()
        st = sd.pop("table")
        ot = od.pop("table")
        return sd == od and np.all(st == ot)

    @property
    def nres(self):
        return len(self.table)

    @property
    def first_MJD(self):
        return self.get_mjds(high_precision=True).min()

    @property
    def last_MJD(self):
        return self.get_mjds(high_precision=True).max()

    def __add__(self, x):
        if type(x) in [int, float]:
            if not x:
                # Adding zero. Do nothing
                return self
        raise NotImplementedError

    def __sub__(self, x):
        if type(x) in [int, float]:
            if not x:
                # Subtracting zero. Do nothing
                return self
        raise NotImplementedError

    def get_mjds(self, high_precision=False):
        """Array of MJDs in the TOAs object

        With high_precision is True
        Return an array of the astropy.times (UTC) of the TOAs

        With high_precision is False
        Return an array of toas in mjd as double precision floats

        WARNING: Depending on the situation, you may get MJDs in a
        different scales (e.g. UTC, TT, or TDB) or even a mixture
        of scales if some TOAs are barycentred and some are not (a
        perfectly valid situation when fitting both Fermi and radio TOAs)
        """
        if high_precision:
            return np.array(self.table["mjd"])
        else:
            return self.table["mjd_float"].quantity

    def get_errors(self):
        """Return a numpy array of the TOA errors in us."""
        return self.table["error"].quantity

    def get_groups(self, gap_limit=None):
        """Flag toas within gap limit (default 2h = 0.0833d) of each other as the same group.

        Groups can be larger than the gap limit - if toas are separated by a gap larger than
        the gap limit, a new group starts and continues until another such gap is found.

        Groups with a two-hour spacing are pre-computed when the TOAs object is constructed,
        and these can rapidly be retrieved from ``self.table`` (which this function will do).

        Parameters
        ----------
        gap_limit : :class:`astropy.units.Quantity`, optional
            The minimum size of gap to create a new group. Defaults to two hours.

        Returns
        -------
        groups : array
            The group number associated to each TOA. Groups are numbered chronologically
            from zero.
        """
        # TODO: make all values Quantity objects for consistency
        if gap_limit is None:
            gap_limit = 2 * u.h
        return _group_by_gaps(self.get_mjds().value, gap_limit.to_value(u.d))
