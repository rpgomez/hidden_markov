#
# Regular cron jobs for the hiddenmarkov package
#
0 4	* * *	root	[ -x /usr/bin/hiddenmarkov_maintenance ] && /usr/bin/hiddenmarkov_maintenance
