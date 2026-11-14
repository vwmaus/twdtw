# Several tests compare the deprecated versions 'f90goto' and 'cpp' of the twdtw
# routine with the default version. The warning that these versions are
# deprecated is given once per session and is tested in test-twdtw_version.R.
# It is marked here as already given, so that it does not appear in other tests.
assign("deprecated_version_warned", TRUE, envir = twdtw:::.twdtw_session)
