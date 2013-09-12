import os
from datetime import datetime

defaultPrefix = os.path.expanduser('~/bin/')

buildTime = datetime.today()
productVersion = 'APBS Trunk '# + buildTime.strftime('%Y%m%d%H%M')

