"""========================================
jnpack     A C++ package manager
Author:    Jian Wang
Email:     wj_hust08@hust.edu.cn
Usage:     python jnpack.py <install|update>
========================================"""

import json
import os
import sys
import glob
import fnmatch

def read_depends(dir):
    path = dir + '/jnpack.json'
    if not os.path.exists(path): return []

    f = open(path)
    content = f.read()
    f.close()

    config = json.loads(content)
    depends = config['depends']
    return depends

def add_packs(ls, n):
    if len(ls) == 0: return

    packs = []
    for pack in ls:
        if not os.path.exists(pack) or n == 0:
            packs.append(pack)
        if not os.path.exists(pack):
            os.system('git clone https://github.com/hust220/%s.git' % (pack,))

    ls = []
    for pack in packs:
        ls.extend(read_depends(pack))

    add_packs(ls, n + 1)

def find_cpps(prefix, dir):
    matches = []
    for root, dirnames, filenames in os.walk(dir):
        for filename in fnmatch.filter(filenames, '*.cpp'):
            matches.append(os.path.join(prefix, root, filename))
    return matches

def update_cmake():
    prefix = os.getcwd()
    ls = []
    for it in glob.glob('*'):
        if it != "CMakeLists.txt":
            ls.extend(find_cpps(prefix, it + '/' + it))
    f = open('CMakeLists.txt', 'w')
    #f.write('add_library(jnpacks STATIC %s)' % (' '.join(ls),))
    f.write('set(jnpacks_sources %s PARENT_SCOPE)' % (' '.join(ls),))
    f.close()

def install():
    depends = read_depends('.')

    if not os.path.exists('jnpacks'):
        os.mkdir('jnpacks')

    os.chdir('jnpacks')
    add_packs(depends, 0)
    update_cmake()
    os.chdir('..')

def update():
    os.chdir('jnpacks')
    for it in os.listdir('.'):
        if it != 'CMakeLists.txt':
            print it
            os.chdir(it)
            os.system('git pull')
            os.chdir('..')
    os.chdir('..')
    os.system('git pull')

def help():
    print __doc__

if __name__ == '__main__':
    if len(sys.argv) == 1:
        help()
    else:
        func = sys.argv[1]
        if func == 'install':
            install()
        elif func == 'update':
            update()
        else:
            help()

