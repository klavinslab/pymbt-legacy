# nupack 3.0 python wrapper. Loosely based on the nupack 2.0 wrapper distributed
# as part of the RBS Calculator:
# (https://github.com/hsalis/Ribosome-Binding-Site-Calculator)
#
# GPL License:
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from subprocess import call
from tempfile import mkdtemp
from shutil import rmtree
from os.path import isdir
from os import environ
from os.path import exists

class nupack:
    """nupack python wrapper class"""
    def __init__(self,sequence_list,material):
        # Set up nupack environment variable
        if 'NUPACKHOME' in environ: 
            self.nupackdir = environ['NUPACKHOME']
        else:
            raise EnvironmentError('NUPACKHOME environment variable must be set') 

        # Store reused information
        if type(sequence_list) is not list:
            if type(sequence_list) is str:
                # if user input a string, make it a list for processing later
                sequence_list = [sequence_list]
            else:
                raise ValueError('Must supply a sequence string or list of them.')
        self.sequence_list = sequence_list
        # Note - the input/output files are not cleaned up by default
        # Use the '_cleanup' method to delete the entire temp dir
        self.outdir = mkdtemp()

        # Set material and check alphabet
        self.material=material
        if material=='rna' or material=='rna1999':
            alphabet = 'AUGCaugc'
        elif material=='dna':
            alphabet = 'ATGCatgc'
        else:
            raise ValueError('material must be \'dna\', \'rna\', or \'rna1999\'.') 

        for i in sequence_list:
            for j,char in enumerate(i):
                if char not in alphabet:
                    raise ValueError('Non-ATGCU character found at index ' + str(j))

        # Keeping track of what commands have been run
        self.complexes_run = False

    def complexes(self,max_complexes,mfe=True,rmdir=True,T=50):
        """
        Runs 'complexes' command on the inputted sequence list.

        'max_complexes' sets the maximum complex size that is calculated (integer).
        'mfe' determines whether mfe calculations are included (boolean, defaults to True).
        'T' sets the temperature.
        """

        # This version can only handle max_complex 2 or below-I think
        self._checkdir() # Recreate temp dir if it was removed

        # Prepare input file
        complexes_input = '\n'.join([str(len(self.sequence_list)), \
                                     '\n'.join(self.sequence_list), \
                                     str(max_complexes)])
        
        f = open(self.outdir+'/nupack.in','w')
        f.write(complexes_input)
        f.close()

        # Prepare arguments
        if mfe==True:
            mfe = ' -mfe '
        else:
            mfe = ''

        # Calculate complexes
        complexes_args = ' -T ' + str(T) + ' -material ' + self.material + ' ' + mfe
        self._run_nupack('complexes',complexes_args)

        # Parse the output
        f_out = open(self.outdir+'/nupack.cx','r+')
        eq_results = f_out.readlines()

        eq_results = [x for x in eq_results if '%' not in x ]
        eq_results = [v.split() for v in eq_results]
       
        # Separate results into complex type and energy
        # Remove rank entry
        for i in eq_results:
            i.pop(0)

        eq_results_cx = [ [ int(float(x.pop(0))) for i,v in enumerate(self.sequence_list) ] for x in eq_results ]
        eq_results_en = [ float(x.pop(0)) for x in eq_results ]
       
        self.complexes_run = max_complexes
        return({'complexes':eq_results_cx, \
                'complex_energy':eq_results_en})

    def concentrations(self,max_complexes,conc=0.5e-6,mfe=True,rmdir=True):
        """
        Runs 'concentrations' command on the inputted sequence list.

        'conc' sets the concentration of each species. Defaults to 0.5e-6 for all.
          Can be either a single value (all species at the same concentraiton) or
          a list of the concentration for each species.
        'max_complexes' sets the maximum complex size that is calculated (integer).
        'mfe' determines whether mfe calculations are included (boolean, defaults to True).
        """
        if self.complexes_run and self.complexes_run is max_complexes:
            pass # complexes has already been run
        else:
            self._checkdir() # Recreate temp dir if it was removed
            # Run complexes in order to get input file
            complexes = self.complexes(max_complexes=max_complexes,mfe=True,rmdir=False)

        # Prepare input file
        if type(conc) == list:
            if len(conc) != len(self.sequence_list):
                raise ValueError("There must be as many 'conc' values as sequences")
            concstring = '\n'.join(conc)
        else:
            concstring = '\n'.join([str(conc) for x in self.sequence_list])

        f = open(self.outdir+'/nupack.con','w')
        f.write(concstring)
        f.close()

        # Run 'complexes'
        concentrations_args = '-sort 3 '
        self._run_nupack('concentrations',concentrations_args)

        # Parse the output of 'complexes'
        f_out = open(self.outdir+'/nupack.eq','r+')
        eq_results = f_out.readlines()
        eq_results = [x for x in eq_results if '%' not in x ]
        eq_results = [v.split() for v in eq_results]
       
        # Format results: complex type, concentration, and energy
        # Remove rank information
        for i in eq_results:
            i.pop(0)

        eq_results_cx = [ [ int(float(x.pop(0))) for i,v in enumerate(self.sequence_list) ] for x in eq_results ]
        eq_results_en = [ x.pop(0) for x in eq_results ]
        eq_results_conc = [ float(x[0]) for x in eq_results ]

        return( {'types':eq_results_cx, \
                'concentration':eq_results_conc, \
                'energy':eq_results_en} )

    def mfe(self,strand,T=50):
        """
        Runs 'mfe' command on the inputted sequence list.

        'T' sets the temperature.
        """

        self._checkdir() # Recreate temp dir if it was removed

        if type(strand) != int:
            raise TypeError('strand value must be integer')
        elif strand > len(self.sequence_list):
            raise ValueError('strand value exceeds the number of sequences')

        sequence = self.sequence_list[strand]
        # Prepare input file
        f = open(self.outdir+'/nupack.in','w')
        f.write(sequence)
        f.close()

        # Run 'mfe'
        mfe_args = ' -T ' + str(T) + ' -material ' + self.material
        self._run_nupack('mfe',mfe_args)

        # Parse the output of 'mfe'
        mfe = float(open(self.outdir+'/nupack.mfe','r+').readlines()[14].strip())

        # Return the mfe of the designated strand
        return(mfe)

    def pairs(self,strand,T=50):
        """
        Runs 'pairs' command on the inputted sequence list.

        'T' sets the temperature.
        """

        self._checkdir() # Recreate temp dir if it was removed

        sequence = self.sequence_list[strand]
        f = open(self.outdir+'/nupack.in','w')
        f.write(sequence)
        f.close()

        # Calculate pair probabilities with 'pairs'
        pairs_args = '-T ' + str(T) + ' -material ' + self.material
        self._run_nupack('pairs',pairs_args)

        # Parse the 'pairs' output
        # Only look at the last n rows - unbound probabilities
        pairs = open(self.outdir+'/nupack.ppairs','r+').readlines()
        pairs = pairs[-len(sequence):]
        pairs = [ x for x in pairs if '%' not in x ]
        #########################
        # Why was this code here?
        #pairs.pop(0) # Remove gap
        #pair_n = pairs.pop(0)
        #########################
        types = [ (int(v.split()[0]),int(v.split()[1])) for v in pairs]
        pp = [ float(v.split()[2]) for v in pairs]

        # Return the pair probabilities
        return({'type':types,'probabilities':pp})

    def _checkdir(self):
        # Create the input data file in a temp dir if it doesn't exist
            if not isdir(self.outdir):
                self.outdir = mkdtemp()

    def _cleanup(self):
        rmtree(self.outdir)

    def _run_nupack(self,command,arguments):
        accepted_commands=['complexes','concentrations','mfe','pairs']
        if command not in accepted_commands:
            raise ValueError('Command must be one of the following:' + accepted_commands)

        env_line = "export NUPACKHOME="+self.nupackdir+' && '
        command_line = 'cd ' + self.outdir + ' && ' + self.nupackdir + '/bin/' + command + ' '
        arguments_line = arguments + ' ' + self.outdir + '/nupack'

        run_command = env_line+command_line+arguments_line
        call(run_command+' > /dev/null',shell=True)

