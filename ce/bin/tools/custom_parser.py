import os
import sys
from ast import literal_eval as leval
from collections import defaultdict
from .environment import log_message

class CustomParser(object):

    current = None

    def __init__(self, cmdline):
        self._raw = cmdline
        self._args = {}

        if CustomParser.current == None:
            CustomParser.current = self

        self.parse_args()

        # return self._args, self._raw[:self._start_index]

    def parse_args(self):

        def replace_dash(value):

            # Top-level pass-in, the value is the nested dict
            if isinstance(value, dict):

                value=dict((key, value) for key, value \
                    in replace_dash(value.iteritems()))

                # This is the final updated dictionary
                yield value

            else:

                # Here, value will either be an iteritems object
                # or a tuple of dict values
                for setting, setting_val in value:

                    setting = setting.replace('-','')

                    if isinstance(setting_val, dict):

                        # Update the nested dictionary by sending an iteritems object
                        setting_val=dict((k,s) for k, s \
                            in replace_dash(setting_val.iteritems()))

                    else:
                        # Just in case, remove dashes from setting-value
                        setting_val = setting_val.replace('-', '')

                        try:
                            leval(setting_val)
                        except:
                            log_message('Could not evaluate custom setting:'
                                ' {}'.format(setting_val))
                        else:
                            setting_val = leval(setting_val)

                        # This returns back to the nested call above
                        yield(setting,setting_val)


                    # This returns back to the top-level pass-in
                    yield setting, setting_val

        args = None

        # Try to find the start of the args
        try:

            if isinstance(self._raw, list):
                start_index = self._raw.index('---all_settings')
                args = self._raw[start_index+1:]

            elif isinstance(self._raw, basestring):
                self._raw = self._raw.split()
                start_index = self._raw.index('---all_settings')
                args = self._raw[start_index+1:]

        except ValueError:

            log_message('Could not parse out custom settings.'
                ' Using default.', depth=2)

            return


        self._start_index = start_index

        args_dict = defaultdict(dict)
        start = 0
        i = start

        g_name = None
        g_setting = None
        g_setting_val = None

        while i < len(args):

            if args[i].startswith('--'):
                g_name = args[i]
                i+=1

                if i >= len(args):
                    break

                while not args[i].startswith('--'):

                    if args[i].startswith('-'):
                        g_setting = args[i]
                        i+=1

                        if i >= len(args):
                            break
                        
                        elif not args[i].startswith('-'):
                            g_setting_val = args[i]
                            args_dict[g_name][g_setting] = g_setting_val

                            i+=1

                            if i>=len(args):
                                break
                    else:
                        i+=1

                        if i>=len(args):
                            break
            else:
                i+=1

        self._args = next(replace_dash(args_dict))

    def retrieve(self, value, default=None):

        if value in self._args:
            return self._args[value]
        else:
            return default

    @classmethod
    def get(cls, value, default=None):
        if cls.current:
            return cls.current.retrieve(value, default)

        else:
            raise RuntimeError('CustomParser has not been initilized')

    @classmethod
    def args(cls, default=None):

        if not cls.current._args:
            return default

        else:
            return cls.current._args

    @classmethod
    def update(cls, genotyper, settings):
        client_args = cls.current.get(genotyper)

        if client_args is None:
            log_message('No client args provided for genotyper: {}'.format(
                genotyper),
                depth=1
            )

            log_message('Using default:')

            for key, value in settings.iteritems():
                log_message('{} -> {}'.format(key, value), depth=2)

            return
        
        log_message('Loading client custom args for genotyper: {}'.format(
            genotyper),
            depth=1
        )

        for key, value in client_args.iteritems():
            settings[key] = value
            log_message('{} -> {}'.format(key, value), depth=2)

if __name__ == '__main__':

    import sys

    test_args = ['---all_settings',
                    '--resistance',
                        '-min', '10.0',
                        '-max', '20',
                        '-test', '100',
                    '--virulence',
                        '-min', '10',
                        '-max', '1000',
                        '-test', '100',
                    '--plasmids',
                        '-min', '10',
                        '-max', 'True',
                        '-test', 'True',
                        ]


    parser = CustomParser(test_args)

    print(parser._args)