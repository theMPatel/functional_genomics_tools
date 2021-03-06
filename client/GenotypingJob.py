#//$autorun;event=PostCreateMainWin

###################################################################
#
# Genotyping Job plugin for all organisms for Bionumerics
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

#imports
import bns
import re
import os
import gzip
import json
import base64
import hashlib
import traceback
import cStringIO
from datetime import datetime
import xml.etree.cElementTree as ET

from collections import defaultdict
from functools import partial, wraps

from wgMLST_Client.SchedulerCommonTools import CeStore
from wgMLST_Client.SchedulerCommonTools.ExecutableJob import ExecutableJob, SingleEntryExecutableJob
from wgMLST_Client.SchedulerCommonTools.CalculationEngineCommunicator import CeCommunicator
from wgMLST_Client.SchedulerCommonTools.CalculationEngineRequest import GetJobResults
from wgMLST_Client.SchedulerCommonTools.ExecJobSettingsDialog import ExecJobSettingsDlg

from wgMLST_Client.Executables.SRSExecutableJob import SRSExecutableJob

from wgMLST_Client.wgMLSTSchema import GetCurrentSchema
from wgMLST_Client.CuratorDbItf import CuratorDbItf
from wgMLST_Client.CommonTools.Settings import StoredSettings

Dlg = bns.Windows.XmlDlgBuilder
MessageBox = bns.Util.Program.MessageBox
ConfirmBox = bns.Util.Program.ConfirmBox

_SALMONELLA = 'Salmonella'
_ESCHERICHIA = 'Escherichia'
_LISTERIA = 'Listeria'
_CAMPY = 'Campylobacter'

organism_abbreviations = {
    'EC': _ESCHERICHIA,
    'SALM': _SALMONELLA,
    'LMO' : _LISTERIA,
    'CAMPY' : _CAMPY
}

default_settings = {

    # ***IMPORTANT***
    # It is incredibly important that the default values
    # here are reflective of the types of input people want:
    # boolean vs float vs integer
    # the types are used explicitly to build specific inputs
    # for the UX
    # 
    _ESCHERICHIA : {
        "plasmids": {
            "percent_identity": 90.0,
            "min_relative_coverage": 90.0,
            "min_merge_overlap": 90.0,
            "search_fragments": False
        },
        "virulence": {
            "percent_identity": 90.0,
            "min_relative_coverage": 60.0,
            "min_merge_overlap": 90.0,
            "search_fragments": True
        },     
        "resistance" : {
            "percent_identity": 90.0,
            "min_relative_coverage": 60.0,
            "min_merge_overlap": 90.0,
            "search_fragments": True
        },
        "ecoli.pathotype" : {
            "percent_identity": 90.0,
            "min_relative_coverage": 60.0,
            "min_merge_overlap": 90.0,
            "search_fragments": True
        },
        "ecoli.serotype" : {
            "percent_identity": 90.0,
            "min_relative_coverage": 60.0,
            "min_merge_overlap": 90.0,
            "search_fragments": False
        },
        "insilicopcr": {
            "max_mismatches": 1,
            "max_nonIUPAC": 2,
            "percent_identity": 90.0,
            "max_length_deviation": 10.0
        },
        # Just use the default settings for now
        "ecoli.stx" :{},
        'fastqc': {
            ''
        }
    },
    _SALMONELLA: {
        "plasmids": {
            "percent_identity": 90.0,
            "min_relative_coverage": 90.0,
            "min_merge_overlap": 90.0,
            "search_fragments": False
        },
        "resistance" : {
            "percent_identity": 90.0,
            "min_relative_coverage": 60.0,
            "min_merge_overlap": 90.0,
            "search_fragments": True
        },
        "salmonella.serotype" : {
            "percent_identity" : 80.0,
            "min_coverage": 70.0,
            "discrimination" : 2.0
        },
        "insilicopcr": {
            "max_mismatches": 3,
            "max_nonIUPAC": 3,
            "percent_identity": 70.0,
            "max_length_deviation": 10.0
        }
    },
    _CAMPY : {
        "plasmids": {
            "percent_identity": 90.0,
            "min_relative_coverage": 90.0,
            "min_merge_overlap": 90.0,
            "search_fragments": False
        },
        "resistance" : {
            "percent_identity": 90.0,
            "min_relative_coverage": 60.0,
            "min_merge_overlap": 90.0,
            "search_fragments": True
        }
    },
    _LISTERIA : {
        "plasmids": {
            "percent_identity": 90.0,
            "min_relative_coverage": 90.0,
            "min_merge_overlap": 90.0,
            "search_fragments": False
        },
        "resistance" : {
            "percent_identity": 90.0,
            "min_relative_coverage": 60.0,
            "min_merge_overlap": 90.0,
            "search_fragments": True
        }
    }
}

needs_validation = {
    'percent_identity',
    'min_coverage',
    'discrimination',
    'min_relative_coverage',
    'min_merge_overlap',
    'max_length_deviation'
}

settings_names = {
    'percent_identity' : 'Percent Identity (%)',
    'discrimination' : 'Discrimination (%)',
    'min_coverage' : "Minimum Coverage (%)",
    'min_relative_coverage' : 'Minimum Relative Coverage (%)',
    'min_merge_overlap': 'Minimum Overlap for Fragment Merge (%)',
    'search_fragments' : 'Search for fragments',
    'max_mismatches' : "Maximum nucleotide mismatches",
    'max_nonIUPAC' : 'Maximum non-IUPAC nucleotides',
    'max_length_deviation': 'Maximum amplicon length deviation (%)'
}

genotyper_names = {
    'plasmids': 'Plasmids',
    'virulence': 'Virulence',
    'resistance': 'Resistance',
    'ecoli.pathotype': 'Pathotype',
    'ecoli.serotype': 'Serotype',
    'salmonella.serotype' : 'Serotype',
    'insilicopcr': 'Insilico PCR',
    'ecoli.stx' : 'Stx subtyper' 
}

# EXPER TYPES
_RESISTANCE_EXPR = 'resistance'
_ANTIBIOTIC_EXPR = 'antibio'
_PLASMIDS_EXPR = 'plasmids'
_VIRULENCE_EXPR = 'virulence'
_STX_EXPR = 'stx'

# RESULTS FILES -> EXPER
_BLAST_RESISTANCE_EXPR = 'resistance.json'
_BLAST_MUTATION_EXPR = 'resistance.point.json'
_BLAST_ANTIBIOTIC_EXPR = 'resistance.antibios.json'

_BLAST_PLASMIDS_EXPR = 'plasmids.json'
_BLAST_VIRULENCE_EXPR = 'virulence.json'
_BLAST_INSILICOPCR_EXPR = 'insilicopcr.json'

_BLAST_ECOLI_PATHOTYPE_EXPR = 'ecoli.pathotype_genotypes.json'
_READS_ECOLI_STX_EXPR = 'ecoli.stxfinder_genotypes.json'
_CONDENSED_ECOLI_STX_EXPR = 'stx_condenser_expr.json'

# FLD NAMES
_SALM_SEROTYPE_FLDNAME = 'Serotype_wgs'
_SALM_ANTIGENIC_FORM_FLDNAME = 'AntigenForm_wgs'
_ECOLI_SEROTYPE_FLDNAME = 'Serotype_wgs'
_ECOLI_TOXIN_FLDNAME = 'Toxin_wgs'
_ECOLI_PATHOTYPE_FLDNAME = 'Pathotype'

# RESULTS FILES -> FIELDS
_CONDENSED_ECOLI_STX_FLD = 'stx_condenser_flds.json'
_READS_SALM_SEROTYPE_FLD = 'salmonella.serotype.json'
_BLAST_ECOLI_PATHOTYPE_FLD = 'ecoli.pathotype_pathotypes.json'
_BLAST_ECOLI_SEROTYPE_FLD = 'ecoli.serotype.json'

results_to_chars = {
    _SALMONELLA: {
        _RESISTANCE_EXPR: [
            _BLAST_RESISTANCE_EXPR,
            _BLAST_MUTATION_EXPR
        ],
        _ANTIBIOTIC_EXPR: [
            _BLAST_ANTIBIOTIC_EXPR
        ],
        _PLASMIDS_EXPR: [
            _BLAST_PLASMIDS_EXPR
        ]
    },
    _ESCHERICHIA: {
        _VIRULENCE_EXPR: [
            _BLAST_VIRULENCE_EXPR,
            _BLAST_ECOLI_PATHOTYPE_EXPR,
            _READS_ECOLI_STX_EXPR,
            _BLAST_INSILICOPCR_EXPR
        ],
        _RESISTANCE_EXPR: [
            _BLAST_RESISTANCE_EXPR,
            _BLAST_MUTATION_EXPR
        ],
        _ANTIBIOTIC_EXPR: [
            _BLAST_ANTIBIOTIC_EXPR
        ],
        _PLASMIDS_EXPR: [
            _BLAST_PLASMIDS_EXPR
        ],
        _STX_EXPR: [
            _CONDENSED_ECOLI_STX_EXPR
        ]
    },
    _CAMPY: {
        _RESISTANCE_EXPR: [
            _BLAST_RESISTANCE_EXPR,
            _BLAST_MUTATION_EXPR
        ],
        _ANTIBIOTIC_EXPR: [
            _BLAST_ANTIBIOTIC_EXPR
        ],
        _PLASMIDS_EXPR: [
            _BLAST_PLASMIDS_EXPR
        ]
    },
    _LISTERIA: {
        _RESISTANCE_EXPR: [
            _BLAST_RESISTANCE_EXPR,
            _BLAST_MUTATION_EXPR
        ],
        _ANTIBIOTIC_EXPR: [
            _BLAST_ANTIBIOTIC_EXPR
        ],
        _PLASMIDS_EXPR: [
            _BLAST_PLASMIDS_EXPR
        ]
    }
}

results_to_fields = {
    
    _SALMONELLA : {
        _SALM_SEROTYPE_FLDNAME : _READS_SALM_SEROTYPE_FLD,
        _SALM_ANTIGENIC_FORM_FLDNAME : _READS_SALM_SEROTYPE_FLD
    },
    _ESCHERICHIA : {
        _ECOLI_PATHOTYPE_FLDNAME : _BLAST_ECOLI_PATHOTYPE_FLD,
        _ECOLI_SEROTYPE_FLDNAME : _BLAST_ECOLI_SEROTYPE_FLD,
        _ECOLI_TOXIN_FLDNAME : _CONDENSED_ECOLI_STX_FLD
    },
    _LISTERIA : {},
    _CAMPY : {}
}

# **** REALLY IMPORTANT ***
# The mapping to genotyper and its specific setting
# is created at runtime by joining the name of the genotyper
# with the setting and using the joined string as the attr id
# when creating a BioNumerics GUI element. It is later parsed and 
# split to get the value. Thus:
# 
# insilicopcr-max_iupac_mismatch as the ID of a GUI element

attr_joiner = '-'

class TabPagesDlg(Dlg.Dialogs):

    def __init__(self, organism):

        Dlg.Dialogs.__init__(self, 'genotyper_settings_dlg')
        self._organism = organism
        self.build_stored_settings()
        self.build_dlg()

    def build_dlg(self):
        # Builds the dialog pages for each of the genotypers
        # and its settings

        def dlg_input(genotyper, setting, default, typ, style, OnChange=None):
            # Here tp (for type) could be anything. Like for example monkey.
            # just think about that

            if style == 'input':
                return Dlg.Input(
                    attr_joiner.join((genotyper, setting)),
                    default=default,
                    tp=typ,
                    remembersettings=True

                )

            elif style == 'radio':

                return Dlg.Radio(
                    attr_joiner.join((genotyper, setting)),
                    ['True', 'False'],
                    30,
                    vertical = False,
                    default = default,
                    remembersettings = True
                )

            elif style == 'check':

                return Dlg.Check(
                    attr_joiner.join((genotyper, setting)),
                    'Activated: ',
                    default = default,
                    remembersettings=True,
                    OnChange=OnChange
                )

        def validate_input(input_ctrl=None, name=''):
            strval = input_ctrl.GetValue()
            try:
                floatval = float(strval)
            except:
                raise RuntimeError('{} input must be a number. '.format(
                    name))

            if floatval < 0.0:
                raise RuntimeError('{} input must be greater than 0.0 '.format(
                    name))

            if floatval > 100.0:
                raise RuntimeError('{} input must be less than 100.0 '.format(
                    name))

        def change_state(args, dlg=None, active_ctrl='', input_ctrls=None):
            # Use the Dialogs class to get the value
            # of the check, which is either 0 or 1
            # as a *string*
            state_value = dlg.GetValue(active_ctrl)

            for ctrl in input_ctrls:
                ctrl.Enabled = bool(int(state_value))

        # Save the inputs for later
        self._inputs = []
        
        tab_dlgs = []

        for genotyper, genotyper_settings in default_settings[self._organism].iteritems():

            grid_here = []
            validations_here = []

            for ux_input_setting, default_setting in genotyper_settings.iteritems():

                if isinstance(default_setting, bool):

                    bn_input_ctrl = dlg_input(
                        genotyper, ux_input_setting, default_setting, None, 'radio')

                elif isinstance(default_setting, int):

                    bn_input_ctrl = dlg_input(
                        genotyper, ux_input_setting, default_setting, 'integer', 'input')

                elif isinstance(default_setting, float):

                    bn_input_ctrl = dlg_input(
                        genotyper, ux_input_setting, default_setting, 'float', 'input')

                grid_here.append([settings_names[ux_input_setting], bn_input_ctrl])

                if ux_input_setting in needs_validation:
                    validations_here.extend([
                        Dlg.Validation(
                            Dlg.ValueOf(bn_input_ctrl).IsNotEmpty(),
                            '{} cannot be empty!'.format(
                                settings_names[ux_input_setting])
                        ),

                        Dlg.ThrowValidate(
                            partial(
                                validate_input,
                                bn_input_ctrl,
                                settings_names[ux_input_setting]
                            ))])

                self._inputs.append(bn_input_ctrl)

            # WARNING: ***NOT IMPLEMENTED***
            # If a user deactivates a certain genotyper
            # we should deactivate the input controls
            # to make sure that it is extra clear
            # that the genotyper is deactivated
            # Otherwise the settings get saved
            # 
            # ***WARNING***
            # We are only protected because the attr active gets set to
            # false so CE ignores that module.

            # state_controller = partial(
            #         change_state,
            #         dlg = self,
            #         active_ctrl='='.join((genotyper, 'activated')),
            #         input_ctrls=tuple(grid_item[1] for grid_item in grid_here)
            #     )


            # state_controller.__name__ = change_state.__name__

            # This check box is so the user can deactivate the genotyper
            active_ctrl = dlg_input(
                genotyper, 'activated', '1', None, 'check')

            # Grid is looking for list of lists
            # Insert this at the front because we want it to show
            # up at the top of the tab page
            grid_here.insert(0, [active_ctrl])

            self._inputs.append(active_ctrl)

            tab_here = Dlg.TabPage(
                grid_here, 
                genotyper_names[genotyper], 
                validations=validations_here
            )

            tab_dlgs.append(tab_here)

        built = Dlg.TabDialog(tab_dlgs, 
            '{} settings'.format(
            self._organism), onOk=self.on_ok)

        dlg = Dlg.Dialog(
            'dlg123',
            '{} settings'.format(self._organism),
            built,
            isresizeable = True
        )

        self.AddDialog(dlg)

    def on_ok(self, args):

        # Get all of the settings after the dialog box has
        # finished
        all_settings = {}

        for input_ctrl in self._inputs:

            str_val = input_ctrl.GetValue()
            attribute = input_ctrl.Id

            if attribute.endswith('activated'):
                val = input_ctrl.Checked
                all_settings[attribute] = bool(val)

                continue

            str_val = input_ctrl.GetValue()

            if attribute.endswith('search_fragments'):
                all_settings[attribute] = str_val == 'True'
                continue

            try:

                if input_ctrl.tp == 'float':
                    actual_val = float(str_val)

                elif input_ctrl.tp == 'integer':
                    actual_val = int(str_val)

                else:
                    raise RuntimeError('Type: {}'.format(input_ctrl.tp))

            except Exception as e:
                genotyper, setting = attribute.split(attr_joiner)

                MessageBox('Error', 'There was an error parsing your'
                    ' custom setting for genotyper:'
                    ' {} and setting: {}.'
                    ' Using default.'.format(
                        genotyper_names[genotyper],
                        settings_names[setting]
                        ) + '\n' + \
                    str(e),
                    '')

                actual_val = default_settings[genotyper][setting]

            all_settings[attribute] = actual_val
        
        bns_settings = StoredSettings(self._organism.upper(), **all_settings)
        bns_settings.Save()

    def build_stored_settings(self):
        all_settings = {}
        bns_settings = StoredSettings(self._organism.upper())
        bns_settings.Load()

        if not len(bns_settings.keys):

            for genotyper, genotyper_settings in \
                default_settings[self._organism].iteritems():

                for parameter_name, default_value in genotyper_settings.iteritems():

                    name = attr_joiner.join([genotyper, parameter_name])
                    all_settings[name] = default_value

                name = attr_joiner.join([genotyper, 'activated'])
                all_settings[name] = True

            bns_settings = StoredSettings(self._organism.upper(), **all_settings)
            bns_settings.Save()

class GenotypingJob(SingleEntryExecutableJob, ExecutableJob):

    _classID = 'Genotyping'
    _displayName = "Genotyping"
    classDescription = "Performing Genotyping"
    
    def __init__(self, entry=None, jobid = '', expertype = None):
        
        if expertype is None:
            expertype = GetCurrentSchema().DeNovoExperTypeName
            
        # specific settings for the executable (required)
        self.acceptanceSettings = StoredSettings(
            'GENOTYPING',
            # **** YOU NEED TO CHANGE THIS!!!! ******
            # **** for the algo you want to run *****
            algorithm='genotyping'
        )
        self.acceptanceSettings.Load()
        
        # generic executable job settings (required)
        GenotypingExecJobSettings = StoredSettings(
            'GENOTYPING_SETTINGS', 
            wallClockTime= 14400, 
            deadlineTime= '',
            requiredMemory = 0, 
            encoding = 'base64', 
            encrypt= False
        )

        GenotypingExecJobSettings.Load()
        
        # create the base class
        SingleEntryExecutableJob.__init__(
            self,
            'Genotyping.bat',
            GenotypingJob.classDescription,
            jobid,
            self.acceptanceSettings,
            GenotypingExecJobSettings,
            expertype,
            entry
        )

        self.define_constants()
        self._org_abbrv = CuratorDbItf.GetCuratorItf().GetOrganismAbbrev()
        
        self._organism = organism_abbreviations.get(self._org_abbrv, None)

        assert self._organism is not None

        self.validate_run()

        self._settings_dlg = TabPagesDlg(self._organism)

        # This will be true when the user runs the genotyping
        # script for the first time
        if not hasattr(self.acceptanceSettings, 'organism'):
            self.acceptanceSettings = StoredSettings(
                'GENOTYPING',
                algorithm='genotyping',
                organism=self._organism
            )

            self.acceptanceSettings.Save()

        elif self.acceptanceSettings.organism != self._organism:
            self.acceptanceSettings = StoredSettings(
                'GENOTYPING',
                algorithm='genotyping',
                organism=self._organism
            )

            self.acceptanceSettings.Save()

    def define_constants(self):
        self.char_mapping = {
                'Not Screened':[0,0],
                'Present': [1,1],
                'Absent': [2,2],
                'Retired': [3,3]
                }

        self.char_present_v = 1.
        self.char_absent_v = 2.
        self.char_retired_v = 3.

        # self.new_char_max = max(x for l in self.char_mapping.itervalues() for x in l)
        # 
        # At first I thought it would be reasonable to just have the max value be the
        # max of the mapping. However, the thought occurred to me that old characters would
        # only get installed with the highest value *at the time*. if we added new mappings
        # the old characters would not be installed with a number that was large enough.
        # 100 should be more than enough for these purposes
        self.new_char_max = 100
    
    def validate_run(self):
        # These are the expertypes that we need

        if self._organism is None:
            raise RuntimeError('You cannot run this job '
                'for the current organism: {}'.format(self._org_abbrv))

        required_chars = results_to_chars[self._organism].keys()
        required_flds = results_to_fields[self._organism].keys()

        missing_chrs = False
        for req in required_chars:
            if not req in bns.Database.Db.ExperimentTypes:
                missing_chrs = True
                break


        missing_flds = False
        for req in required_flds:
            if req not in bns.Database.Db.Fields:
                missing_flds = True
                break

        if missing_flds or missing_chrs:
            answer = ConfirmBox('Missing database fields/experiments for job:'
                ' Genotyping. Would you like to install?')

            if answer == 1:
                self.install_dependencies(required_chars, required_flds)
                return True

            else:
                return False

    def install_dependencies(self, chars, flds):

        chrs_to_add = []

        for chrname in chars:
            if chrname not in bns.Database.Db.ExperimentTypes:
                chrs_to_add.append(chrname)

        for expr in chrs_to_add:
            et = bns.Database.Db.ExperimentTypes.Create(expr,
                'CHR', displayName=expr)

            et.Settings = '<Settings><Settings><Settings OpenSet="1"'\
            ' AbsetIsZero="1"/></Settings></Settings>'

            chrType = bns.Characters.CharType(expr)

            for name, mapRange in self.char_mapping.iteritems():
                chrType.MapAdd(name,mapRange[0], mapRange[1])

        flds_to_add = []
        for fldname in flds:
            if fldname not in bns.Database.Db.Fields:
                flds_to_add.append(fldname)


        for fld in flds_to_add:
            bns.Database.Db.Fields.Add(fld)

        return True

    def _GetSubmitClientVersion(self):
        """Submit job XML version - default 1
        BLAST allele finder version 2 : guarantees that the k-mer size is present in the allowed list
        of the BLAST algorithm.
        """
        return '1'

    def show_settings_dlg(self, args):
        self._settings_dlg.Show()

    #for display in dialog
    def ShowOnDialog(self, row, validations, withCellText=True):

        self._settings_button = Dlg.Button( 'action', 'Settings', contextCall=\
            Dlg.ContextCall(self.show_settings_dlg), buttonID='settings_button' )

        row.append([self._organism, self._settings_button])
    
        return True
    
    #this is called just before the dialog opens
    def InitDlgValues(self, dlg):
        return
    
    # this is called after the user pressed OK in the dlg box
    def ReadFromDlg(self, dlg):
        return

    # actually creates the Dlg box
    def GetSettingsDlg(self):
        return ExecJobSettingsDlg('WgMLSTClient' + self._classID + 'SettingsDlg', self)
    
    # This is where you add whatever you want to add to the cmdline
    def GetExtraCmdLine(self):
        # note: params that are given on the constructor to the base class, are automagically on the cmdline

        cmdline = ''
        reads = self.get_readfiles()

        if reads:

            for read in reads:
                cmdline += '--query-reads {} '.format(read)
        # cmdline_dict = self.ce_genotyper_settings()

        # cmdline = '---all_settings '

        # for genotyper, settings in cmdline_dict.iteritems():
        #     str_here = ' '.join('-{arg} {param}'.format(arg=key, param=value) for \
        #         key, value in settings.iteritems())
        #     cmdline += '--{0} {1} '.format(genotyper, str_here)

        localArgsDict = {}
        ceComm = CeCommunicator.GetCurrent()
        if ceComm.GetCommunicationVersion() == '1':
            localArgsDict['--query'] = '"[RESULTSDIR]\\denovo.fasta.gz"'

        return cmdline + ' '.join('--{0} {1}'.format(key, value) for key, value in localArgsDict.iteritems())

    def ce_genotyper_settings(self):
        # Note that if this object (self) changes at runtime, 
        # your updates will not be kept which is why I have to do the below.
        # s = StoredSettings('GENOTYPING')
        # s.Load()

        # If the above were not the case, I could do this and expect it
        # to work:
        # 
        # settings_key = self.acceptanceSettings.organism

        # s.organism is the organism

        # settings_key = s.organism.upper()
        settings_key = self._organism.upper()

        bns_settings = StoredSettings(settings_key)
        bns_settings.Load()

        cmdline_dict = defaultdict(dict)

        for attr in bns_settings.keys:

            genotyper, setting = attr.split(attr_joiner)

            if isinstance(getattr(bns_settings, attr), float):
                cmdline_dict[genotyper][setting] = str(getattr(bns_settings, attr)/100.)
            else:
                cmdline_dict[genotyper][setting] = str(getattr(bns_settings, attr))

        self.job_dependencies(self._organism, cmdline_dict)

        return cmdline_dict
    
    #this is where you send along files
    def GetFiles(self):

        # Return as tuples?
        files_to_return = [(
            'genotyper_settings.json',
            base64.b64encode(json.dumps(self.ce_genotyper_settings())),
            False
        )]

        ceComm = CeCommunicator.GetCurrent()
        if ceComm.GetCommunicationVersion() != '1':
            return files_to_return
        
        if not bns.Database.Experiment(self.Entry, self.ExperType).IsPresent():
            raise RuntimeError("No assembled sequence present for entry '{0}' and experiment type '{1}'. BLAST job canceled.".format(self.Entry.Key, self.ExperType.Name))

        seq = bns.Sequences.SequenceData(self.Entry.Key, self.ExperType.Name).LoadSequence()
        if not seq:
            raise RuntimeError("Internal error: no denovo sequence for entry {0}.".format(self.Entry.Key))
        
        flName = 'denovo.fasta.gz'
        fgz = cStringIO.StringIO()
        with gzip.GzipFile(flName, 'wb', fileobj=fgz) as fl:
            self.ExportSeqToFasta(seq, fl)
        
        seqCompr = base64.b64encode(fgz.getvalue())
        encrypted = False

        files_to_return.append((flName, seqCompr, encrypted))

        return files_to_return

    def ExportSeqToFasta(self, seq, fl):
        # do *not* change the header format ('denovo_NNN')
        # it is used verbatim to parse the contig number in later
        # processing, e.g. in the wgMLST QA window
        for contigIdx, contig in enumerate(seq.split('|')):
            c = re.sub('[^ACGTUWSMKRYBDHV]','N', contig.upper()) # replace everythings that's not a iupac code with N
            fl.write('>denovo_{0}\n{1}\n'.format(contigIdx, c))
        
    def GetDataList(self, calcComm=None):
        dataList = super(GenotypingJob, self).GetDataList(calcComm)
        
        ceComm = CeCommunicator.GetCurrent()
        if ceComm.GetCommunicationVersion() != '1':
            dataList.append({
                        'DataID': self.GetDeNovoFileId(calcComm),
                        'Type': 'Upload',
                        'Prefix': '--query '
                    })

        return dataList
    
    def GetDeNovoFileId(self, calcComm):
        seq = bns.Sequences.SequenceData(self.Entry.Key, self.ExperType.Name).LoadSequence()
        if not seq:
            raise RuntimeError("Internal error: no denovo sequence for entry {0}.".format(self.Entry.Key))
        
        m =  hashlib.md5()
        m.update(seq)
        denovoHash = m.hexdigest()
        
        settingHash = 'CE_DENOVOHASH'
        settingFileId = 'CE_DENOVOFILEID'
        expAttach = bns.Database.ExperAttachments(self.Entry.ID, self.ExperType.ID)
        denovoFileId = expAttach.get(settingFileId)
        
        ceComm = CeCommunicator.GetCurrent()
        if denovoFileId and ceComm.CheckUploadFileId(denovoFileId) and denovoHash == expAttach.get(settingHash):
            return denovoFileId
        
        flName = 'denovo.fasta.gz'
        flPath = os.path.join(bns.Database.Db.Info.TempDir, flName)
        with gzip.GzipFile(flPath, 'wb') as fl:
            self.ExportSeqToFasta(seq, fl)
        
        denovoFileId = ceComm.UploadFile(flPath, calcComm)
        expAttach[settingFileId] = denovoFileId
        expAttach[settingHash] = denovoHash
        
        return denovoFileId

    def get_readfiles(self):

        reads = bns.Database.Experiment(self.Entry.Key, 'wgs').LoadOrCreate()
        network_files = reads.GetLinks()
        to_return = []

        splitters = {
            'PulseNetSneakerTest' : '\\sneakerTest',
            'PulsenetSneakerDev' : '\\sneakerDev',
            'PulsenetSneakerProd' : '\\sneakerTest'
        }

        if len(network_files) != 2:
            return []

        for i in xrange(len(network_files)):
            for splitter in splitters:

                if splitter in network_files[i]:
                    to_keep = network_files[i].split(splitter)[1]
                    to_keep = list(filter(None, to_keep.split(os.sep)))
                    to_return.append(os.path.join(splitters[splitter], *to_keep))
                    break

        if len(to_return) != 2:
            return []

        else:
            return to_return
    
    #add a follow-up job
    def AddJob(self, executableJob):
        pass

    def job_dependencies(self, organism, cmdline_dict):
        # Use this function to modify the cmdline_dictionary
        # with anything else you want to plop onto the CE

        def salmonella():

            if 'salmonella.serotype' not in cmdline_dict:
                return

            cmdline_dict['salmonella.serotype']['ani_value'] = 'None'

            # For the future in case that we want to accept
            # precomputed ANI values
            # 
            # # Make sure the ANI field exists
            # if 'ANI' not in bns.Database.Db.Fields:
            #     cmdline_dict['salmonella.serotype']['ani_value'] = None

            # else:
            #     # Add the value to the dictionary
            #     ani_value = bns.Database.EntryField(self.Entry, 'ANI').Content
            #     cmdline_dict['salmonella.serotype']['ani_value'] = ani_value

        def escherichia():
            return

        def campy():
            return

        def listeria():
            return

        funcs = {
            'Salmonella' : salmonella,
            'Escherichia' : escherichia,
            'Campylobacter': campy,
            'Listeria' : listeria
        }
        # Call the appropriate value
        funcs[organism]()
    
    # Process job results
    def Process(self, calcComm, xml=None):
        #fetch the results
        jobDisplay = "result of {0} of entry {1}".format(GenotypingJob.GetExecutableDisplayName(), self.Entry.Key)
        requester = GetJobResults(self)
        bns.Util.Program.AddActionLog("Start processing {0}".format(jobDisplay))
        results = requester.DoRequest(calcComm)
        bns.Util.Program.AddActionLog("Got result from server for {0}".format(jobDisplay))
        
        try:
            self._Process(results)
        except:
            MessageBox('', traceback.format_exc(), '')
            raise
        
        bns.Util.Program.AddActionLog("Done processing  {0}".format(jobDisplay))
        
        return True
    
    # Implementation of the processing
    def _Process(self, results):

        if not len( results ):
            return True

        for key in results:
            try:
                results[key] = json.loads(base64.b64decode(results[key]))
            except:
                raise RuntimeError('Could not parse genotyping results for entry: {}'.format(
                    self.Entry.Key))

        def field_processer(organism, to_process):

            if not isinstance(to_process, dict):
                raise RuntimeError('Inappropriate results '
                    'information of type: {}'.format(type(to_process)))

            def salmonella():
                for fld, information in to_process.iteritems():

                    if not 'results' in information or not \
                        len(information['results']):
                        continue

                    to_field = ''

                    if fld == _SALM_SEROTYPE_FLDNAME:
                        to_field = information['results'].get('serotype', '')

                    elif fld == _SALM_ANTIGENIC_FORM_FLDNAME:
                        to_field = information['results'].get('formula', '')

                    self.Entry.Field(fld).Content = to_field

                bns.Database.Db.Fields.Save()

            def escherichia():
                for fld, information in to_process.iteritems():
                    
                    if not 'results' in information or not \
                        len(information['results']):
                        continue

                    to_field = ''

                    if fld == _ECOLI_PATHOTYPE_FLDNAME:

                        to_field = information['results'][0]
                    
                    elif fld == _ECOLI_SEROTYPE_FLDNAME:

                        format_str = '{otype}:{htype}'

                        otypes = set(information['results']['O'])
                        htypes = set(information['results']['H'])

                        if not len(otypes) and not len(htypes):
                            to_field = ''
                        
                        else:
                            to_field = format_str.format(
                                otype = '/'.join(otypes),
                                htype = '/'.join(htypes)
                            )


                    elif fld == _ECOLI_TOXIN_FLDNAME:
                        to_field = '; '.join(information['results'])

                    self.Entry.Field(fld).Content = to_field

                bns.Database.Db.Fields.Save()

            def listeria():
                pass

            def campy():
                pass

            organisms = {
                _SALMONELLA: salmonella,
                _ESCHERICHIA: escherichia,
                _LISTERIA : listeria,
                _CAMPY : campy
            }

            if not organisms.get(organism, False):
                raise NotImplementedError(
                    'Processing for organism:' 
                    '{} not implemented'.format(
                    organism))

            else:
                organisms[organism]()

        # testing = os.path.expanduser('~\\Desktop\\testing.json')

        # with open(testing, 'w') as f:
        #     json.dump(results, f)

        chars_to_results = results_to_chars.get(self._organism, {})
        flds_to_results = results_to_fields.get(self._organism, {})

        for expr, file_names in chars_to_results.iteritems():

            bn_expr = bns.Characters.CharSetType(expr)
            found_chars = { x:self.char_absent_v for x in xrange(bn_expr.GetCount()) }
            retired_chars = set(map(bn_expr.GetChar, found_chars.iterkeys()))
            new_chars = {}

            for file_name in file_names:

                results_here = results.get(file_name, None)

                if results_here is None:
                    continue

                if 'results' not in results_here:
                    continue

                chars_here = results_here['results']

                for character, presence in chars_here.iteritems():

                    charNr = bn_expr.FindChar(character)
                    retired_chars.discard(character)

                    if charNr > 0:
                        found_chars[charNr] = self.char_present_v if presence else self.char_absent_v

                    else:
                        new_chars[character] = presence

            if len(new_chars):
                bn_expr.NoSave()

                for new in new_chars:
                    charNr = bn_expr.AddChar(new, int(self.new_char_max))
                    found_chars[charNr] = self.char_present_v if new_chars[new] else self.char_absent_v

                bn_expr.Save()

            if len(retired_chars):

                indices = map(bn_expr.FindChar, retired_chars)

                for i in indices:
                    found_chars[i] = self.char_retired_v

            vector = [found_chars[i] for i in xrange(bn_expr.GetCount())]

            this_entry_expr = bns.Database.Experiment(self.Entry, expr).LoadOrCreate()
            this_entry_expr.SetData(vector)
            this_entry_expr.Save()

        to_process = defaultdict(list)
        for fld, file_name in flds_to_results.iteritems():
            if not file_name in results:
                continue

            to_process[fld] = results[file_name]

        if to_process:
            field_processer(self._organism, to_process)

        # Save the results as a json dump
        attachment_obj = self.get_attachment('Genotyping results')
        attachment_obj.Content = json.dumps(results)

        # Build all the pretty tables that are needed:
        tables_needed = {
            'ResFinder results' : 'resistance.json'
        }

        for table_name, results_name in tables_needed.iteritems():

            if results_name in results:
                self.build_pretty_table(table_name, results[results_name])

        return True
    
    def build_pretty_table(self, table_name, results):

        if not 'extra' in results:
            return

        attachment_obj = self.get_attachment(name=table_name)
        
        table_header = [
            'locus',
            'contig_id',
            'identity',
            'coverage',
            'query_start',
            'query_stop',
            'full_match'
        ]

        # Make the header
        table = '\t'.join(table_header) + '\n'

        extra_results = results['extra']

        for result in extra_results:

            line = []
            if not 'hits' in result or not \
                len(result['hits']):
                continue

            if len(result['hits']) == 1:
                for heading in table_header:

                    if heading in result:
                        line.append(result[heading])

                    else:
                        line.append(result['hits'][0][heading])

                for i in xrange(len(line)):

                    if not isinstance(line[i], basestring):

                        if isinstance(line[i], float):
                            line[i] = str(round(line[i], 4))

                        else:
                            line[i] = str(line[i])

                table += '\t'.join(line) + '\n'

            else:

                for hit in result['hits']:

                    new_line = []

                    for heading in table_header:

                        if heading in result:
                            new_line.append(result[heading])

                        else:
                            new_line.append(hit[heading])
                    
                    for i in xrange(len(new_line)):

                        if not isinstance(new_line[i], basestring):

                            if isinstance(new_line[i], float):
                                new_line[i] = str(round(new_line[i], 4))

                            else:
                                new_line[i] = str(new_line[i])

                    line.append(new_line)

                for l in line:

                    table += '\t'.join(l) + '\n'

        attachment_obj.Content = table

    def get_attachment(self, name=''):

        if not name:
            return

        to_return = None

        for attachment in self.Entry.GetAttachmentList():
            if attachment['ClassID'] == 'TEXT' and \
                attachment['Name'] == name:

                to_return = attachment['Attach']
                break

        if to_return is None:
            to_return = self.Entry.AddAttachmentFlatText(
                '', True, name=name)


        return to_return

    #used to check whether this job can be run on a particular entry
    def IsStartExpPresent(self):
        return bns.Database.Experiment(self.Entry, self.ExperType).IsPresent()

        
#job registration
from wgMLST_Client.SchedulerCommonTools.ExecutableJob import ExecutableJob
ExecutableJob.RegisterExecutable(GenotypingJob)

#testing
if __name__ == '__main__' and __bnsdebug__:
    import bns
    #from wgMLST_Client import initcontext
    print ExecutableJob.GetExecutables()
    
    r = GenotypingJob()
    bns.Stop()
    
    from wgMLST_Client.WgmlstLaunchSingleEntryJobsDlg import WgmlstLaunchSingleEntryJobsDlg
    dlg = WgmlstLaunchSingleEntryJobsDlg('wgMLSTClientLaunchJobsDlg', bns.Database.Db.Selection, 1, {}, '')
    if dlg.Show():
        print len(dlg.toDoJobs)
        if len(dlg.toDoJobs):
            print dlg.toDoJobs[0].GetCommandLine()

