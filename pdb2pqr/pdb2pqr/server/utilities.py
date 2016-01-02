""" Helper functions for web-based use of PDB2PQR """

def startLogFile(jobName, fileName, logInput):
    with open('%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobName, fileName), 'w') as f:
        f.write(logInput)

def appendToLogFile(jobName, fileName, logInput):
    with open('%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobName, fileName), 'a') as f:
        f.write(logInput)

def resetLogFile(jobName, fileName):
    """
    For clearing out old log files if needed.
    Used mainly for removing apbs_end_time if apbs is rerun.
    """
    filename = '%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobName, fileName)
    try:
        os.remove(filename)
    except EnvironmentError:
        pass

def getTrackingScriptString(jobid=None):
    """
    For injecting tracking script into a web page.

    jobid -> current jobid. Adds "jobid" custom variable to events and page views on this page.
    """
    customVarString = ""

    if jobid is not None:
        customVarString = "_gaq.push(['_setCustomVar',1,'jobid','{jobid}',3]);".format(jobid=str(jobid))

    #If you look closely you'll see escaped { and }.
    string = """<script type="text/javascript">

  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-11026338-3']);
  _gaq.push(['_setDomainName', 'none']);
  _gaq.push(['_setAllowLinker', true]);
  {customVar}
  _gaq.push(['_trackPageview']);

  (function() {{
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  }})();

</script>""".format(customVar=customVarString)
    return string

def getEventTrackingString(category, action, label, value=None):
    valueString = ', {value}'.format(value=value) if value is not None else ""
    eventString = '_gaq.push(["_trackEvent", "{category}", "{action}", "{label}"{valuestr}]);\n'
    return eventString.format(category=str(category), action=str(action), label=str(label), valuestr=valueString)

def createHTMLTypeMap(protein, definition, outfilename):
    """
        Create an HTML typemap file at the desired location. If a
        type cannot be found for an atom a blank is listed.

        Parameters
            protein:  A Protein Object
            definition: The definition objects.
            outfilename:  The name of the file to write (string)
    """
    from .forcefield import Forcefield
    from aconf import STYLESHEET

    # Cache the initial atom numbers
    numcache = {}
    for atom in protein.getAtoms():
        numcache[atom] = atom.serial
    protein.reSerialize()

    amberff = Forcefield("amber", definition, None)
    charmmff = Forcefield("charmm", definition, None)

    file = open(outfilename, "w")
    file.write("<HTML>\n")
    file.write("<HEAD>\n")
    file.write("<TITLE>PQR Typemap (beta)</TITLE>\n")
    file.write("<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">\n" % STYLESHEET)
    file.write("</HEAD>\n")
    file.write("<BODY>\n")
    file.write("<H3>This is a developmental page including the atom type for the atoms in the PQR file.</H3><P>\n")
    file.write("<TABLE CELLSPACING=2 CELLPADDING=2 BORDER=1>\n")
    file.write("<tr><th>Atom Number</th><th>Atom Name</th><th>Residue Name</th><th>Chain ID</th><th>AMBER Atom Type</th><th>CHARMM Atom Type</th></tr>\n")

    for atom in protein.getAtoms():
        if isinstance(atom.residue, (Amino, WAT, Nucleic)):
            resname = atom.residue.ffname
        else:
            resname = atom.residue.name

        ambergroup = amberff.getGroup(resname, atom.name)
        charmmgroup  = charmmff.getGroup(resname, atom.name)


        file.write("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n" % (atom.serial, atom.name, resname, atom.chain_id, ambergroup, charmmgroup))


    file.write("</table>\n")
    file.write("</BODY></HTML>\n")
    file.close()

    # Return the original numbers back
    for atom in protein.getAtoms():
        atom.serial = numcache[atom]

    del numcache
    del amberff
    del charmmff
