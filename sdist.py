"""
Specialization of the Python source distributor to copy stuff to standard location. This
only needs to work for me.
"""

import os, shutil, datetime
from distutils.command.sdist import sdist as _sdist

class sdist(_sdist):

    def make_distribution(self):
        # standard stuff
        _sdist.make_distribution(self)

        # extra: copies tar.gz file over to software web pages, updates link
        # and the index.html file
        if 'WEB_PATH' in os.environ:
            dest_dir = os.path.join(os.environ['WEB_PATH'], 'software')
            if os.path.isdir(dest_dir):
                # name of tar file with version number
                full_name = self.distribution.get_fullname() + '.tar.gz'

                # tar file with version number and directory added
                base_name = os.path.join(self.dist_dir, full_name)

                if os.path.isfile(base_name):
                    # name without version number
                    oname = self.distribution.get_name()
                    name  = oname + '.tar.gz'

                    # name with version number in web page directory
                    dest_name = os.path.join(dest_dir, full_name)
                    print 'copying',base_name,'->',dest_name
                    shutil.copy(base_name, dest_name)

                    # soft link
                    link_name = os.path.join(dest_dir, name)
                    print 'soft linking',dest_name,'->',link_name
                    if os.path.islink(link_name):
                        os.remove(link_name)
                    if os.path.exists(link_name):
                        print link_name,'exists but should not; please fix.'
                        exit(1)
                    os.symlink(dest_name, link_name)


                    # html file
                    index_file = os.path.join(dest_dir, 'index.html')
                    if os.path.isfile(index_file):
                        print 'updating',index_file
                        fin = open(index_file)
                        version = self.distribution.get_version()
                        lines = []
                        now = datetime.datetime.now()
                        for line in fin:
                            if line.find(oname + ' version') > -1:
                                line = '<!-- ' + oname + ' version --><td>' + version + '</td>\n'
                            if line.find(oname + ' date') > -1:
                                line = '<!-- ' + oname + ' date --><td>' + now.ctime() + '</td>\n'
                            lines.append(line)
                        fin.close()

                        fout = open(index_file, 'w')
                        for line in lines:
                            fout.write(line)
                        fout.close()

                        
