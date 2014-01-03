Import('env')

def check_custom(context, program, message, exit_on_fail=True):
    context.Message( message )
    ret = context.TryCompile(program, '.c')
    if ret == 0:
        if exit_on_fail:
            print "test failed!\n"
            print "Error: cannot run the following test program:"
            print program
            print "Please check your compiler flags."
            Exit(1)
        return None
    ret = context.TryRun(program, '.c')
    context.Result(ret[0])
    return ret

def checkEPSILON(context, test_type):
    program = """
#include <stdio.h>
int main()
{{
    {0} epsilon = 1.0;
    while(1.0 + epsilon / 2.0 != 1.0)
    {{
        epsilon /= 2.0;
    }}
    printf( "%e", epsilon );
    return 0;
}}
""".format(test_type)

    return check_custom(context, program, 'getting size of ' + test_type + '... ')[1]

def checkEMBED(context):
    program = """
#define EMBED(rctag)

static const char* rctag;
static void* use_rcsid=(0 ? &use_rcsid : (void**)&rcsid);
EMBED(rcsid)
"""

    result = check_custom(context, program, 'checking if preprocessor can embed macros...', False)
    return result is not None

conf = Configure(env, custom_tests = {'checkEPSILON' : checkEPSILON, 'checkEMBED' : checkEMBED})

float_tag = '@FLOAT_EPSILON@'
double_tag = '@DOUBLE_EPSILON@'
defined = ''
undefined = '//'

config_dict = {}
config_dict[float_tag] = conf.checkEPSILON('float')
print "float epsilon is: " + config_dict[float_tag] 
config_dict[double_tag] = conf.checkEPSILON('double')
print "double epsilon is: " + config_dict[double_tag] 

can_embed = conf.checkEMBED()
print "can embed macros: " + str(can_embed)
config_dict['@USE_MACRO_EMBED@'] = defined if can_embed else undefined

config_dict['@HAVE_TIME@'] = defined if conf.CheckFunc('time') else undefined
config_dict['@HAVE_RAND@'] = defined if conf.CheckFunc('rand') else undefined
config_dict['@HAVE_SRAND@'] = defined if conf.CheckFunc('srand') else undefined

config_dict['@USE_READLINE@'] = defined if conf.CheckLib('readline') else undefined
config_dict['@USE_ZLIB@'] = defined if conf.CheckLib('z') else undefined

if not conf.CheckLib('m'):
    print "libm is required..."
    Exit(1)

config_dict['@USE_TINKER@'] = undefined
config_dict['@USE_FETK@'] = undefined
config_dict['@USE_MPI@'] = undefined
config_dict['@USE_DEBUG@'] = undefined
config_dict['@USE_VERBOSE_DEBUG@'] = undefined
config_dict['@USE_VAPBSQUIET@'] = undefined

config_dict['@PACKAGE_STRING@'] = '1.4'

env.Substfile('src/apbscfg.h.in', SUBST_DICT = config_dict)

#maloc configuration stuff.
config_dict['@ACCEPT_USES_UINT@'] = undefined
config_dict['@ACCEPT_USES_ULONG@'] = undefined
config_dict['@HAVE_CYGWIN@'] = undefined
config_dict['@HAVE_ARPA_INET_H@'] = defined if conf.CheckCHeader('arpa/inet.h')  else undefined
config_dict['@HAVE_FCNTL_H@'] = defined if conf.CheckCHeader('fcntl.h')  else undefined
config_dict['@HAVE_GETCWD@'] = defined if conf.CheckFunc('getcwd')  else undefined
config_dict['@HAVE_INTTYPES_H@'] = defined if conf.CheckCHeader('inttypes.h')  else undefined
config_dict['@HAVE_MEMORY_H@'] = defined if conf.CheckCHeader('memory.h')  else undefined

if False:
    config_dict['@HAVE_MPI_H@'] = defined if conf.CheckLibWithHeader('inttypes.h')  else undefined
else:
    config_dict['@HAVE_MPI_H@'] = undefined
    
config_dict['@HAVE_NETDB_H@'] = defined if conf.CheckCHeader('netdb.h')  else undefined

config_dict['@HAVE_O_NONBLOCK@'] = defined if conf.CheckDeclaration('O_NONBLOCK', '#include <fcntl.h>\n')  else undefined

config_dict['@HAVE_READLINE_HISTORY_H@'] = defined if conf.CheckCHeader('readline/history.h')  else undefined
config_dict['@HAVE_READLINE_READLINE_H@'] = defined if conf.CheckCHeader('readline/readline.h')  else undefined

if conf.CheckCHeader('rpc/xdr.h'):
    config_dict['@HAVE_RPC_XDR_H@'] = config_dict['@HAVE_XDR@'] = defined
else:
    config_dict['@HAVE_RPC_XDR_H@'] = config_dict['@HAVE_XDR@'] = undefined
config_dict['@HAVE_NETINET_IN_H@'] = defined if conf.CheckCHeader('netinet/in.h')  else undefined
config_dict['@HAVE_STDINT_H@'] = defined if conf.CheckCHeader('stdint.h')  else undefined
config_dict['@HAVE_STDLIB_H@'] = defined if conf.CheckCHeader('stdlib.h')  else undefined
config_dict['@HAVE_STRINGS_H@'] = defined if conf.CheckCHeader('strings.h')  else undefined
config_dict['@HAVE_STRING_H@'] = defined if conf.CheckCHeader('string.h')  else undefined
config_dict['@HAVE_SYS_SOCKET_H@'] = defined if conf.CheckCHeader('sys/socket.h')  else undefined
config_dict['@HAVE_SYS_STAT_H@'] = defined if conf.CheckCHeader('sys/stat.h')  else undefined

config_dict['@HAVE_SYS_TIMES_H@'] = defined if conf.CheckCHeader('sys/times.h')  else undefined
config_dict['@HAVE_SYS_TIME_H@'] = defined if conf.CheckCHeader('sys/time.h')  else undefined
config_dict['@HAVE_SYS_TYPES_H@'] = defined if conf.CheckCHeader('sys/types.h')  else undefined
config_dict['@HAVE_SYS_UN_H@'] = defined if conf.CheckCHeader('sys/un.h')  else undefined
config_dict['@HAVE_SYS_WAIT_H@'] = defined if conf.CheckCHeader('sys/wait.h')  else undefined
config_dict['@HAVE_UNISTD_H@'] = defined if conf.CheckCHeader('unistd.h') else undefined

config_dict['@MALOC_PACKAGE_NAME@'] = 'maloc'
config_dict['@MALOC_PACKAGE_VERSION@'] = '1.0.1'
config_dict['@MALOC_PACKAGE_STRING@'] = config_dict['@MALOC_PACKAGE_NAME@'] + ' ' + config_dict['@MALOC_PACKAGE_VERSION@']
config_dict['@MALOC_PACKAGE_TARNAME@'] = config_dict['@MALOC_PACKAGE_STRING@']

config_dict['@STAT_MACROS_BROKEN@'] = undefined

std_headers = ('stdlib.h','stdarg.h', 'string.h', 'float.h')
has_std_headers = all((conf.CheckCHeader(h) for h in std_headers))

config_dict['@STDC_HEADERS@'] = defined if has_std_headers else undefined

config_dict['@MODE_T_INT@'] = defined if not conf.CheckType('mode_t', '#include <sys/types.h>\n') else undefined
config_dict['@PID_T_INT@'] = defined if not conf.CheckType('pid_t', '#include <sys/types.h>\n') else undefined
config_dict['@SIZE_T_UINT@'] = defined if not conf.CheckType('size_t', '#include <sys/types.h>\n') else undefined

config_dict['@CRT_SECURE_NO_DEPRECATE@'] = undefined
config_dict['@CRT_NONSTDC_NO_DEPRECATE@'] = undefined

env.Substfile('src/maloc/src/maloccf.h.in', SUBST_DICT = config_dict)
