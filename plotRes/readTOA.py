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

from pint.observatory import Observatory, get_observatory
from pint.pulsar_mjd import Time

__all__ = [
    "TOAs",
    "get_TOAs",
    "format_toa_line",
]

toa_commands = (
    "DITHER",
    "EFAC",
    "EMAX",
    "EMAP",
    "EMIN",
    "EQUAD",
    "FMAX",
    "FMIN",
    "INCLUDE",
    "INFO",
    "JUMP",
    "MODE",
    "NOSKIP",
    "PHA1",
    "PHA2",
    "PHASE",
    "SEARCH",
    "SIGMA",
    "SIM",
    "SKIP",
    "TIME",
    "TRACK",
    "ZAWGT",
    "FORMAT",
    "END",
)


def get_TOAs(timfile):
    """Load and prepare TOAs for PINT use.
    """
    t = TOAs(timfile)
    return t

def _toa_format(line, fmt="Unknown"):
    """Determine the type of a TOA line.

    Identifies a TOA line as one of the following types:
    Comment, Command, Blank, Tempo2, Princeton, ITOA, Parkes, Unknown.
    """
    if re.match(r"[0-9a-z@] ", line):
        return "Princeton"
    elif (
        line.startswith("C ")
        or line.startswith("c ")
        or line[0] == "#"
        or line.startswith("CC ")
    ):
        return "Comment"
    elif line.upper().startswith(toa_commands):
        return "Command"
    elif re.match(r"^\s+$", line):
        return "Blank"
    elif re.match(r"^ ", line) and len(line) > 41 and line[41] == ".":
        return "Parkes"
    elif len(line) > 80 or fmt == "Tempo2":
        return "Tempo2"
    elif re.match(r"\S\S", line) and len(line) > 14 and line[14] == ".":
        # FIXME: This needs to be better
        return "ITOA"
    else:
        return "Unknown"


def _parse_TOA_line(line, fmt="Unknown"):
    """Parse a one-line ASCII time-of-arrival.

    Return an MJD tuple and a dictionary of other TOA information.
    The format can be one of: Comment, Command, Blank, Tempo2,
    Princeton, ITOA, Parkes, or Unknown.
    """
    MJD = None
    fmt = _toa_format(line, fmt)
    d = dict(format=fmt)
    if fmt == "Princeton":
        # Princeton format
        # ----------------
        # columns  item
        # 1-1     Observatory (one-character code) '@' is barycenter
        # 2-2     must be blank
        # 16-24   Observing frequency (MHz)
        # 25-44   TOA (decimal point must be in column 30 or column 31)
        # 45-53   TOA uncertainty (microseconds)
        # 69-78   DM correction (pc cm^-3)
        d["obs"] = get_observatory(line[0].upper()).name
        d["freq"] = float(line[15:24])
        d["error"] = float(line[44:53])
        ii, ff = line[24:44].split(".")
        MJD = (int(ii), float("0." + ff))
        try:
            d["ddm"] = float(line[68:78])
        except ValueError:
            d["ddm"] = 0.0
    elif fmt == "Tempo2":
        # This could use more error catching...
        fields = line.split()
        d["name"] = fields[0]
        d["freq"] = float(fields[1])
        if "." in fields[2]:
            ii, ff = fields[2].split(".")
            MJD = (int(ii), float("0." + ff))
        else:
            MJD = (int(fields[2]), 0.0)
        d["error"] = float(fields[3])
        d["obs"] = get_observatory(fields[4].upper()).name
        # All the rest should be flags
        flags = fields[5:]
        for i in range(0, len(flags), 2):
            k, v = flags[i].lstrip("-"), flags[i + 1]
            if k in ["error", "freq", "scale", "MJD", "flags", "obs", "name"]:
                raise ValueError(f"TOA flag ({k}) will overwrite TOA parameter!")
            try:  # Convert what we can to floats and ints
                d[k] = int(v)
            except ValueError:
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
    elif fmt == "Command":
        d[fmt] = line.split()
    elif fmt == "Parkes":
        """
        columns     item
        1-1         Must be blank
        26-34       Observing Frequency (MHz)
        35-55       TOA (decimal point must be in column 42)
        56-63       Phase offset (fraction of P0, added to TOA)
        64-71       TOA uncertainty
        80-80       Observatory (1 character)
        """
        d["name"] = line[1:25]
        d["freq"] = float(line[25:34])
        ii = line[34:41]
        ff = line[42:55]
        MJD = (int(ii), float("0." + ff))
        phaseoffset = float(line[55:62])
        if phaseoffset != 0:
            raise ValueError(
                "Cannot interpret Parkes format with phaseoffset=%f yet" % phaseoffset
            )
        d["error"] = float(line[63:71])
        d["obs"] = get_observatory(line[79].upper()).name
    elif fmt == "ITOA":
        raise RuntimeError("TOA format '%s' not implemented yet" % fmt)
    return MJD, d


def format_toa_line(
    toatime,
    toaerr,
    freq,
    obs,
    dm=0.0 * u.pc / u.cm ** 3,
    name="unk",
    flags={},
    format="Princeton",
):
    """Format TOA line for writing

    Parameters
    ----------
    toatime
        Time object containing TOA arrival time
    toaerr
        TOA error as a Quantity with units
    freq
        Frequency as a Quantity with units (NB: value of np.inf is allowed)
    obs
        Observatory object
    dm
        DM for the TOA as a Quantity with units (not printed if 0.0 pc/cm^3)
    name
        Name to embed in TOA line (conventionally the data file name)
    format
        (Princeton | Tempo2)
    flags
        Any Tempo2 flags to append to the TOA line

    Returns
    -------
    out : str
        Formatted TOA line

    Note
    ----
    This implementation does not undo things like ``TIME`` statements; when used
    by :func:`pint.toa.TOAs.write_TOA_file` these commands are not emitted either.

    Princeton format::

        columns  item
        1-1     Observatory (one-character code) '@' is barycenter
        2-2     must be blank
        16-24   Observing frequency (MHz)
        25-44   TOA (decimal point must be in column 30 or column 31)
        45-53   TOA uncertainty (microseconds)
        69-78   DM correction (pc cm^-3)

    Tempo2 format:

        - First line of file should be "``FORMAT 1``"
        - TOA format is ``name freq sat satErr siteID <flags>``
    """
    if format.upper() in ("TEMPO2", "1"):
        toa_str = Time(toatime, format="pulsar_mjd_string", scale=obs.timescale)
        # In Tempo2 format, freq=0.0 means infinite frequency
        if freq == np.inf * u.MHz:
            freq = 0.0 * u.MHz
        flagstring = ""
        if dm != 0.0 * u.pc / u.cm ** 3:
            flagstring += "-dm {0:%.5f}".format(dm.to(u.pc / u.cm ** 3).value)
        # Here I need to append any actual flags
        for flag in flags.keys():
            v = flags[flag]
            # Since toas file do not have values with unit in the flags,
            # here we are taking the units out
            if flag in ["clkcorr"]:
                continue
            if hasattr(v, "unit"):
                v = v.value
            flag = str(flag)
            if flag.startswith("-"):
                flagstring += " %s %s" % (flag, v)
            else:
                flagstring += " -%s %s" % (flag, v)
        # Now set observatory code. Use obs.name unless overridden by tempo2_code
        try:
            obscode = obs.tempo2_code
        except AttributeError:
            obscode = obs.name
        out = "%s %f %s %.3f %s %s\n" % (
            name,
            freq.to(u.MHz).value,
            toa_str,
            toaerr.to(u.us).value,
            obscode,
            flagstring,
        )
    elif format.upper() in ("PRINCETON", "TEMPO"):
        # This should probably use obs.timescale instead of this hack
        if obs.tempo_code == "@":
            toa_str = str(Time(toatime, format="pulsar_mjd_string", scale="tdb"))
        else:
            toa_str = str(Time(toatime, format="pulsar_mjd_string", scale="utc"))
        # The Princeton format can only deal with MJDs that have up to 20
        # digits, so truncate if longer.
        if len(toa_str) > 20:
            toa_str = toa_str[:20]
        # In TEMPO/Princeton format, freq=0.0 means infinite frequency
        if freq == np.inf * u.MHz:
            freq = 0.0 * u.MHz
        if obs.tempo_code is None:
            raise ValueError(
                "Observatory {} does not have 1-character tempo_code!".format(obs.name)
            )
        if dm != 0.0 * u.pc / u.cm ** 3:
            out = obs.tempo_code + " %13s%9.3f%20s%9.2f                %9.4f\n" % (
                name,
                freq.to(u.MHz).value,
                toa_str,
                toaerr.to(u.us).value,
                dm.to(u.pc / u.cm ** 3).value,
            )
        else:
            out = obs.tempo_code + " %13s%9.3f%20s%9.2f\n" % (
                name,
                freq.to(u.MHz).value,
                toa_str,
                toaerr.to(u.us).value,
            )
    else:
        raise ValueError("Unknown TOA format ({0})".format(format))

    return out


def read_toa_file(filename, process_includes=True, cdict=None):
    """Read TOAs from the given filename into a list.
    """
    if isinstance(filename, str):
        with open(filename, "r") as f:
            return read_toa_file(f, process_includes=process_includes, cdict=cdict)
    else:
        f = filename

    ntoas = 0
    toas = []
    commands = []
    if cdict is None:
        cdict = {
            "EFAC": 1.0,
            "EQUAD": 0.0 * u.us,
            "EMIN": 0.0 * u.us,
            "EMAX": np.inf * u.us,
            "FMIN": 0.0 * u.MHz,
            "FMAX": np.inf * u.MHz,
            "INFO": None,
            "SKIP": False,
            "TIME": 0.0,
            "PHASE": 0,
            "PHA1": None,
            "PHA2": None,
            "MODE": 1,
            "JUMP": [False, 0],
            "FORMAT": "Unknown",
            "END": False,
        }
        top = True
    else:
        top = False

    for line in f.readlines():
        MJD, d = _parse_TOA_line(line, fmt=cdict["FORMAT"])
        if d["format"] == "Command":
            cmd = d["Command"][0].upper()
            commands.append((d["Command"], ntoas))
            if cmd == "SKIP":
                cdict[cmd] = True
                continue
            elif cmd == "NOSKIP":
                cdict["SKIP"] = False
                continue
            elif cmd == "END":
                cdict[cmd] = True
                break
            elif cmd in ("TIME", "PHASE"):
                cdict[cmd] += float(d["Command"][1])
            elif cmd in ("EMIN", "EMAX", "EQUAD"):
                cdict[cmd] = float(d["Command"][1]) * u.us
            elif cmd in ("FMIN", "FMAX", "EQUAD"):
                cdict[cmd] = float(d["Command"][1]) * u.MHz
            elif cmd in ("EFAC", "PHA1", "PHA2"):
                cdict[cmd] = float(d["Command"][1])
                if cmd in ("PHA1", "PHA2", "TIME", "PHASE"):
                    d[cmd] = d["Command"][1]
            elif cmd == "INFO":
                cdict[cmd] = d["Command"][1]
                d[cmd] = d["Command"][1]
            elif cmd == "FORMAT":
                if d["Command"][1] == "1":
                    cdict[cmd] = "Tempo2"
            elif cmd == "JUMP":
                if cdict[cmd][0]:
                    cdict[cmd][0] = False
                    cdict[cmd][1] += 1
                else:
                    cdict[cmd][0] = True
            elif cmd == "INCLUDE" and process_includes:
                # Save FORMAT in a tmp
                fmt = cdict["FORMAT"]
                cdict["FORMAT"] = "Unknown"
                new_toas, new_commands = read_toa_file(d["Command"][1], cdict=cdict)
                toas.extend(new_toas)
                commands.extend(new_commands)
                # re-set FORMAT
                cdict["FORMAT"] = fmt
            else:
                continue
        if cdict["SKIP"] or d["format"] in ("Blank", "Unknown", "Comment", "Command"):
            continue
        elif cdict["END"]:
            if top:
                break
        else:
            newtoa = TOA(MJD, **d)
            if (
                (cdict["EMIN"] > newtoa.error)
                or (cdict["EMAX"] < newtoa.error)
                or (cdict["FMIN"] > newtoa.freq)
                or (cdict["FMAX"] < newtoa.freq)
            ):
                continue
            else:
                newtoa.error *= cdict["EFAC"]
                newtoa.error = np.hypot(newtoa.error, cdict["EQUAD"])
                if cdict["INFO"]:
                    newtoa.flags["info"] = cdict["INFO"]
                if cdict["JUMP"][0]:
                    newtoa.flags["jump"] = cdict["JUMP"][1]
                if cdict["PHASE"] != 0:
                    newtoa.flags["phase"] = cdict["PHASE"]
                if cdict["TIME"] != 0.0:
                    newtoa.flags["to"] = cdict["TIME"]
                toas.append(newtoa)
                ntoas += 1

    return toas, commands


def build_table(toas, filename=None):
    mjds, mjd_floats, errors, freqs, obss, flags = zip(
        *[
            (
                t.mjd,
                t.mjd.mjd,
                t.error.to_value(u.us),
                t.freq.to_value(u.MHz),
                t.obs,
                t.flags,
            )
            for t in toas
        ]
    )
    return table.Table(
        [
            np.arange(len(mjds)),
            table.Column(mjds),
            np.array(mjd_floats) * u.d,
            np.array(errors) * u.us,
            np.array(freqs) * u.MHz,
            np.array(obss),
            np.array(flags),
        ],
        names=(
            "index",
            "mjd",
            "mjd_float",
            "error",
            "freq",
            "obs",
            "flags",
        ),
        meta={"filename": filename},
    ).group_by("obs")


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

class TOA:

    def __init__(
        self,
        MJD,
        error=0.0,
        obs="Barycenter",
        freq=float("inf"),
        scale=None,
        flags=None,
        **kwargs,
    ):
        site = get_observatory(obs)
        # If MJD is already a Time, just use it. Note that this will ignore
        # the 'scale' argument to the TOA() constructor!
        if isinstance(MJD, time.Time):
            if scale is not None:
                raise ValueError("scale argument is ignored when Time is provided")
            t = MJD
        else:
            try:
                arg1, arg2 = MJD
            except TypeError:
                arg1, arg2 = MJD, None
            if scale is None:
                scale = site.timescale
            # First build a time without a location
            # Note that when scale is UTC, must use pulsar_mjd format!
            if scale.lower() == "utc":
                fmt = "pulsar_mjd"
            else:
                fmt = "mjd"
            t = time.Time(arg1, arg2, scale=scale, format=fmt, precision=9)

        # Now assign the site location to the Time, for use in the TDB conversion
        # Time objects are immutable so you must make a new one to add the location!
        # Use the intial time to look up the observatory location
        # (needed for moving observatories)
        # The location is an EarthLocation in the ITRF (ECEF, WGS84) frame
        try:
            loc = site.earth_location_itrf(time=t)
        except Exception:
            # Just add informmation and re-raise
            raise
        # Then construct the full time, with observatory location set
        self.mjd = time.Time(t, location=loc, precision=9)

        if hasattr(error, "unit"):
            try:
                self.error = error.to(u.microsecond)
            except u.UnitConversionError:
                raise u.UnitConversionError(
                    "Uncertainty for TOA with incompatible unit {0}".format(error)
                )
        else:
            self.error = error * u.microsecond
        self.obs = site.name
        if hasattr(freq, "unit"):
            try:
                self.freq = freq.to(u.MHz)
            except u.UnitConversionError:
                raise u.UnitConversionError(
                    "Frequency for TOA with incompatible unit {0}".format(freq)
                )
        else:
            self.freq = freq * u.MHz
        if self.freq == 0.0 * u.MHz:
            self.freq = np.inf * u.MHz
        if flags is None:
            self.flags = kwargs
        else:
            self.flags = flags
            if kwargs:
                raise TypeError(
                    f"TOA constructor does not accept keyword arguments {kwargs} when flags are specified."
                )

    def __str__(self):
        s = (
            self.mjd.mjd_string
            + f": {self.error.value:6.3f} {self.error.unit} error at '{self.obs}' at {self.freq.value:.4f} {self.freq.unit}"
        )
        if self.flags:
            s += " " + str(self.flags)
        return s

    def as_line(self, format="Tempo2", name=None, dm=0 * u.pc / u.cm ** 3):
        """Format TOA as a line for a ``.tim`` file."""
        if name is None:
            name = self.name
        return format_toa_line(
            mjd=self.mjd,
            error=self.error,
            freq=self.freq,
            obs=self.obs,
            dm=dm,
            name=name,
            format=format,
            flags=self.flags,
        )


class TOAs:

    def __init__(self, toafile=None, toalist=None):
        # First, just make an empty container
        self.commands = []
        self.filename = None

        toalist, self.commands = read_toa_file(toafile)

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
    def ntoas(self):
        return len(self.table)

    @property
    def observatories(self):
        return set(self.get_obss())

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

    def get_freqs(self):
        """Return a numpy array of the observing frequencies in MHz for the TOAs"""
        return self.table["freq"].quantity

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

    def get_obss(self):
        """Return a numpy array of the observatories for each TOA."""
        return self.table["obs"]

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


    def write_TOA_file(self, filename, name="unk", format="tempo2"):
        """Write this object to a ``.tim`` file.

        This function writes the contents of this object to a (single) ``.tim``
        file. If ``TEMPO2`` format is used, this file is able to represent the
        contents of this object to nanosecond level. No ``TIME`` or ``EFAC``
        commands are emitted.

        Parameters
        ----------
        filename : str or file-like
            File name to write to; can be an open file object
        name : str
            Value to put in the "name" field of tempo2 files, if a "-name" flag is
            not available.
        format : str
            Format specifier for file ('TEMPO' or 'Princeton') or ('Tempo2' or '1');
            note that not all features may be supported in 'TEMPO' mode.
        """
        try:
            # FIXME: file must be closed even if an exception occurs!
            # Answer is to use a with statement and call the function recursively
            outf = open(filename, "w")
            handle = False
        except TypeError:
            outf = filename
            handle = True

        if format.upper() in ("TEMPO2", "1"):
            outf.write("FORMAT 1\n")

        # Add pulse numbers to flags temporarily if there is a pulse number column
        # FIXME: everywhere else the pulse number column is called pulse_number not pn

        for (toatime, toaerr, freq, obs) in zip(
            self.table["mjd"],
            self.table["error"].quantity,
            self.table["freq"].quantity,
            self.table["obs"],
        ):
            obs_obj = Observatory.get(obs)

            toatime_out = toatime
            out_str = format_toa_line(
                toatime_out,
                toaerr,
                freq,
                obs_obj,
                format=format,
            )
                #name=flags.pop("name", name),
            outf.write(out_str)

        if not handle:
            outf.close()
