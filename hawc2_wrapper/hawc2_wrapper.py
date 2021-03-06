""" Wrapper for HAWC2, HAWC2S, and HAWC2aero  """

import platform
import os
import zipfile
import glob
import shutil
import time
try:
    import subprocess32 as subprocess
    _with_timeout = True
except:
    import subprocess
    _with_timeout = False
    print 'Warning: reinstall wrapper or download subprocess32 package!'
import warnings


class HAWC2Wrapper(object):
    """
    Wrapper for HAWC2, HAWC2S, and HAWC2aero.
    It only executes the codes and checks the log files. It does not read any
    result file.

    Parameters
    ----------
    hawc2bin: str
        Name of the executable to run. The executable can be an HAWC2,
        HAWCStab2, and HAWC2aero executable.
    wine_cmd: str
        Name of the wine application to run HAWC2 on a non-windows machine.
    case_id: str
        Name of the htc file to execute.
    data_directory: str
        Name of data directory.
    res_directory: str
        Name of results directory.
    log_directory: str
        Name of logfile directory.
    copyback_results: bool
        Flag to copy results and model to a desired folder.
    copyback_results_dir: str
        Name of folder where to copy the results and model
    dry_run: bool
        Flag to skip the execution.

    Returns
    -------
        The class does not return or create anything.

    Example
    -------
    >>> from hawc2_wrapper import HAWC2Wrapper
    >>> executer = HAWC2Wrapper()
    >>> executer.hawc2bin = 'HAWC2s.exe'
    >>> executer.case_id = 'new_file_h2s'
    >>> executer.execute()

    """

    def __init__(self, **kwargs):
        super(HAWC2Wrapper, self).__init__()

        self.hawc2bin = 'hawc2MB.exe'
        self.wine_cmd = 'wine'
       
        self.case_id = 'hawc2_case'

        self.data_directory = 'data'
        self.res_directory = 'res'
        self.log_directory = ''

        self.copyback_results = True
        self.copyback_results_dir = 'res_copy'

        self.verbose = True

        self.dry_run = False
        self.timeout = 1000

        self.basedir = os.getcwd()

        for k, w in kwargs.iteritems():
            try:
                setattr(self, k, w)
            except:
                pass

    def compute(self):

        if self.verbose: print 'executing %s for case:  %s.' % (self.hawc2bin, self.case_id)

        tt = time.time()

        exec_str = []

        # check platform and add wine to command list if we're on OS X or Linux
        _platform = platform.platform()
        if 'Linux' in _platform or 'Darwin' in _platform:
            exec_str.append(self.wine_cmd)

        exec_str.append(self.hawc2bin)
        exec_str.append(self.case_id+'.htc')

        if not self.dry_run:
            try:
                if _with_timeout:
                    proc = subprocess.check_output(exec_str,
                                                   timeout=self.timeout)
                else:
                    proc = subprocess.check_output(exec_str)
                self.success = True
                if self.verbose: print self.hawc2bin, 'output:'
                if self.verbose: print proc
            except:
                self.success = False
                print self.hawc2bin, ' crashed for case %s' % self.case_id
        else:
            print self.hawc2bin + ' dry run...'
            self.success = True

        if 'hawc2s' in self.hawc2bin.lower():
            self.check_log(['Error', 'iterations exceeded'])

        elif 'hawc2mb' in self.hawc2bin.lower():
            try:
                self.check_log(['Error'])
            except:
                print 'no log file found for case %s ...' % self.case_id

        if self.copyback_results:
            self.copy_results()

        if self.verbose: print 'HAWC2Wrapper simulation time for case %s: %f' % (self.case_id,
                                                                time.time()-tt)

    def copy_results(self):
        """
        copy results files into the folder copyback_results_dir
        """
        results_dir = os.path.join(self.basedir, self.copyback_results_dir)
        try:
            os.mkdir(results_dir)
        except:
            pass

        # grab all files in the folder with name starting with self.case_id
        files = glob.glob(self.case_id + '*.*')
        files.append('error.out')

        # copy files in folder
        for filename in files:
            try:
                shutil.copy(filename, results_dir)
            except:
                if self.verbose: print 'failed copying back file "%s" into %s' %\
                                     (filename, results_dir)

        # copy data directory in folder
        try:
            shutil.rmtree(os.path.join(results_dir, self.data_directory),
                          ignore_errors=True)
            shutil.copytree(self.data_directory,
                            os.path.join(results_dir, self.data_directory))
        except:
            print 'failed copying back data directory for' +\
                                 ' case %s' % self.case_id

        # copy results directory in folder
        if 'hawc2mb' in self.hawc2bin.lower():
            try:
                shutil.rmtree(os.path.join(results_dir, self.res_directory),
                              ignore_errors=True)
                shutil.copytree(self.res_directory,
                                os.path.join(results_dir, self.res_directory))
            except:
                print 'failed copying back res directory for' +\
                                     ' case %s' % self.case_id

    def check_log(self, messages):
        """
        prints lines of the log file containing specific messages
        """
        logfile = os.path.join(self.log_directory, self.case_id + '.log')
        with open(logfile, 'r') as fid:
            for line in fid.readlines():
                for m in messages:
                    if (m in line):
                        warnings.warn(line)

    def zipdir(self, path):
        """
        zip files in a directory
        """
        zip = zipfile.ZipFile(path + '.zip', 'w')

        for root, dirs, files in os.walk(path):
            for file in files:
                zip.write(os.path.join(root, file))

        zip.close()

    def copy_dirs(self, directory, patterns, overwrite=False):
        """
        copy files from execution directory back to simulation root
        """

        if isinstance(patterns, basestring):
            patterns = [patterns]

        for pattern in patterns:
            for src_path in sorted(glob.glob(pattern)):
                dst_path = os.path.join(directory, pattern)
                if overwrite:
                    try:
                        shutil.rmtree(dst_path)
                    except:
                        pass
                print src_path
                try:
                    shutil.copytree(src_path, dst_path)
                except:
                    raise RuntimeError('Copy failed - directory probably \
                                        exists. use overwrite = True')

if __name__ == '__main__':
    pass
