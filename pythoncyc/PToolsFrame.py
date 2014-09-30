"""
Copyright (c) 2014, SRI International

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
----------------------------------------------------------------------

See class PFrame below for most of the documentation about this module.

"""

import sys
from PTools import PToolsError, PythonCycError
import pythoncyc.config as config
if 'IPython' in sys.modules:
    from IPython.display import display, HTML

def convertLispIdtoPythonId(s):
    """
     Convert string s such that it can be a valid Python identifier.
    """
    return s.replace('-','_').replace('.','_').replace('?','_p').replace('|','').lower()


class PFrame():
    """
    PFrame is a class to represent Pathway Tools' frames. A PFrame can
    represent a class frame (e.g., Reactions) as well as an instance frame
    (e.g., RXN-9000). 

    The required parameters to create a PFrame are the frame id (a string) and a PGDB
    object. For example, assuming that 'meta' is bound to a PGDB object, the following
    create a PFrame to represent the reaction RXN-9000,

        PFrame('RXN-9000', meta)

    By default, an instance PFrame (not a class PFrame) is created.
    To create a PGDB object, use class PGDB, or call method pythoncyc.so.

    By default, the data of the frame is not requested from the server, that is,
    a PFrame object is created containing only the frame id and the PGDB.
    All its slots and data can be fetched from Pathway Tools by
    specifying the keyword argument getFrameData=True. For example, the following
    create a frame for reaction RXN-9000 and retrieve all its slots and data,

        PFrame('RXN-9000', meta, getFrameData=True)

    For creating a class, the isClass keyword parameter must say so,

        PFrame('Reactions', meta, isClass=True)

    When creating a class, if getFrameData=True is specified, the class slots and its
    data are fetched and all the instances of the class are also created.

    In all cases, accessing one slot of a PFrame object will bring in all slots data
    for that PFrame, if the data for that PFrame has not been retrieved.
    For example, for a reaction, the slots left, right, substrates, and so on,
    are populated if one of the slots of its PFrame object is accessed.

    The attribute values of the PFrames cannot be modified, that is, slots are read only.
    On the other hand, slots of Pathway Tools' objects can be modified using methods
    put_slot_value and put_slot_values. See class PGDB for these methods.

    IPython
    -------
    
    A few methods were written to display frames using HTML in IPython.
    Naturally, this functionality is only defined when using IPython as a
    Python interpreter.
    
    """

    def __init__(self, frameid, pgdb, getFrameData=False, isClass = False):
        """
        Creation of a PFrame is done lazily when getFrameData is False but is retrieved
        now when getFrameData is True.
        If getFrameData is False, the PFrame is created assuming that the frameid exists in
        the PGDB in which case an access to a slot will trigger the transfer of the
        whole frame.
        """
        # Always store the frameid surrounded by vertical bars.
        self.__dict__['frameid'] =  frameid if (frameid.startswith('|') and frameid.endswith('|')) else '|'+frameid+'|'
        self._isclass = isClass
        if isClass:
           self.__dict__['instances'] = []
        # The complete frame (the slot values) is not created by default. 
        self._gotframe = False
        # This is the PGDB object, not the PGDB name.
        self.__dict__['pgdb'] = pgdb
        # Add this frame on the list of current frames for this pgdb.
        pgdb.__dict__['_frames'][convertLispIdtoPythonId(frameid)] = self
        # TBD: add the frame on the list of instances for the corresponding class.
        if getFrameData:
            self.get_frame_data()
        return None

    # The following four definitions are for the pickle (or cPickle) module.

    def __getinitargs__(self):
        return (self.frameid, self.pgdb)
    
    def __getstate__(self):
        # Use a compact representation of a PFrame.
        return self.vectorize_dict()

    def __setstate__(self, values):
        # TBD: Use schema of this PFrame to recreate dictionary.
        self.__dict__ = dict

    def vectorize_dict(self):
        self.__dict__.values()

    # End of definitions for pickle.

    def __getslice__(self, i, j, stride = None):
        if config._debug:
            print 'PFrame __getslice__ ', i, j, stride
        return self.instances[i:j:stride]

    def __getattr__(self, attr):
        # Accessing a slot of the frame using attribute syntax (e.g. r.left)
        if config._debug:
            print 'PFrame __getattr__ ', attr

        if attr in self.__dict__:
            return self.__dict__[attr]
        else:
            attrId = convertLispIdtoPythonId(attr)
            if attrId in self.__dict__:
                return self.__dict__[attrId]
                
        if self._gotframe:
                return None
                # raise PythonCycError('No slot with name %s exists for frame %s.' % (attr, self.frameid))
        else: 
            # Get the frame object from Pathway Tools.
            self.get_frame_data()
            if not (attrId in self.__dict__):
                return None
                # raise PythonCycError('No slot with name %s exists for frame %s.' % (attr, self.frameid))
            else:
                return self.__dict__[attr]

    def __getitem__(self,attr):
        if config._debug:
           print "PFrame __getitem__ ", attr
        # The slice case is for attr = slice(i,j,s)
        if (isinstance(attr,int) or isinstance(attr, slice)):
           if 'instances' in self.__dict__:
              return self.instances[attr]
           else:
              raise PythonCycError('Indexing cannot be applied because this is a PFrame which is not a vector and this PFrame has no instances attribute.')
        
        if attr in self.__dict__:
            return self.__dict__[attr]
        else:
            attrId = convertLispIdtoPythonId(attr)
            if attrId in self.__dict__:
                return self.__dict__[attrId]
            else:
                return self.__getattr__(attr)

    def __dir__(self):
        return (dir(self.__class__) + self.__dict__.keys())

    def get_frame_data(self):
        """
        Retrieve the frame data from Pathway Tools, that is, all slots and their values
        for this PFrame are retrieved and stored locally. For a class, the instances
        are not retrieved by this method. Instead, use method get_class_data applied to a PGDB
        object.

        Return
           the self PFrame, modified with the new slots and data.
        """
        # FrameObject is a dictionary of slot names and values.
        frameObject = self.pgdb.sendPgdbFnCall('get-frame-object', self.frameid)
        if not frameObject:
            raise PythonCycError("Could not retrieve frame "+self.frameid+" from organism (orgid) "+self.pgdb._orgid)
        else:
            self._gotframe = True
            # Modify slot names to allow Python's syntax (e.g., '_' instead of '-').
            for slot in frameObject:
                self.__dict__[convertLispIdtoPythonId(slot)] = frameObject[slot]
        return self
    
    def __setattr__(self, attr, val):
        if not attr.startswith('_'):
            raise PythonCycError("Attributes of PFrame objects cannot be modified "+str((self.frameid,attr,val))+". Use methods put_slot_value or put_slot_values using a PGDB object instead.")
        self.__dict__[attr] = val
        return None
    
    def __setitem__(self, attr, val):
        if not attr.startswith('_'):
            raise PythonCycError("Attributes of PFrame objects cannot be modified "+str((self.frameid,attr,val))+". Use methods put_slot_value or put_slot_values using a PGDB object instead.")
        self.__dict__[attr] = val
        return None
    
    def __str__(self):
        if self._isclass:
            return '<PFrame class ' + self.frameid+' currently with '+str(len(self.instances))+' instances ('+self.pgdb._orgid+')>'
        else:
            return '<PFrame instance ' + self.frameid+' ('+self.pgdb._orgid+')>'
        
    def __repr__(self):
        if self._isclass:
           return self.__str__()
        else:
           return str(self.__dict__)

    if 'IPython' in sys.modules:
        def _ipython_display_(self):
          table = ''
          if not self._isclass:
              for attr in sorted(self.__dict__):
                  if not (attr.startswith('_')):
                     table = table+"<tr><td>"+str(attr)+"</td><td>"+str(self.__dict__[attr])+"</td></tr>"
          else:
              table = table+'<tr><td>Class '+self.frameid+' has '+str(self._nb_pframes())+' instances</td></td></tr>'
          if table == '':
             table = '<tr><td>' + self.__str__() + '</td></tr>'
          table = '<table>'+ table + '</table>'
          display(HTML(table))

    def _nb_pframes(self):
        return len(self.instances)

    def __eq__(self, frame):
        # If the given frame is a string, it is assumed to be a frameid.
        if isinstance(frame, str):
            return self.frameid == frame
        # Otherwise it could be a Frame.
        elif isinstance(frame, PFrame): 
            # TBD: should we consider the name of the pgdb?
            return (self.pgdb == frame.pgdb) and (self.frameid == frame.frameid)
        # Otherwise it cannot be equal to this object.
        else: return False

    def __ne__(self, frame):
        return not self.__eq__(frame)


