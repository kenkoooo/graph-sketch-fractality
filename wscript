# -*- python -*-
import fnmatch
import os

APPNAME= 'agl'
VERSION= '0.0.1'

top = '.'
out = 'bin'

def options(opt):
  opt.load('compiler_cxx')
  opt.load('unittest_gtest')
  opt.add_option('--build_debug', action='store_true', default=False, help='debug build')
  opt.add_option('--build_profile', action='store_true', default=False, help='debug build')

def configure(conf):
  conf.load('compiler_cxx')
  conf.load('unittest_gtest')
  conf.env.CXXFLAGS += ['-Wall', '-g', '-std=c++0x', '-pthread']
  conf.env.LINKFLAGS += ['-pthread']

  if conf.options.build_debug:
    conf.env.CXXFLAGS += ['-O0', '-D_GLIBCXX_DEBUG']
  elif conf.options.build_profile:
    conf.env.CXXFLAGS += ['-O3', '-fno-inline']
    conf.env.LINKFLAGS += ['-Wl,--no-as-needed', '-lprofiler', '-Wl,--as-needed']
  else:
    conf.env.CXXFLAGS += ['-O3']

def build(bld):
  bld.recurse('3rd_party')

  cc_file_main = []
  cc_file_test = []
  cc_file_stlib = []
  target_dirs = ['src', 'playground']
  if bld.cmd != 'build':
    target_dirs.append('tutorial')

  for src_dirname in target_dirs:
    for root, dirnames, filenames in os.walk(src_dirname):
      for filename in fnmatch.filter(filenames, '*.cc'):
        if filename == 'test.cc':
          print '\033[91m\033[1mWarning: You should avoid the filename \'test.cc\'.'
          print 'It may conflict with gtest and cause unexpected behavior.\033[0m'
        filepath = os.path.join(root, filename)
        if filename.endswith('_main.cc'):
          cc_file_main.append(filepath)
        elif filename.endswith('_test.cc'):
          cc_file_test.append(filepath)
        else:
          cc_file_stlib.append(filepath)

  bld.stlib(
    name     = 'agl',
    source   = cc_file_stlib,
    target   = 'agl',
    includes = ['src', 'playground', '3rd_party'])

  for cc in cc_file_main:
    n = os.path.basename(cc).replace('_main.cc', '')
    bld.program(
      name     = n,
      source   = cc,
      target   = n,
      use      = ['agl', 'gflags'],
      includes = ['src', 'playground', '3rd_party'])

  bld.program(
    name     = 'test',
    source   = cc_file_test,
    target   = 'test',
    use      = ['agl', 'gflags', 'gtest'],
    includes = ['src', 'playground', '3rd_party'])

  bld.program(
    name     = 'testt',
    features = 'testt',
    source   = cc_file_test,
    target   = 'testt',
    use      = ['agl', 'gflags', 'gtest'],
    includes = ['src', 'playground', '3rd_party'])

from waflib.Build import BuildContext
class build_full(BuildContext):
  cmd = 'build-full'
