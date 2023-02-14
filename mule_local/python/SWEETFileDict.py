import struct
import numpy as np



class SWEETFileDict:
    """
    Write binary file with data which can be read again into SWEET.
    
    The file format is as follows:
    
    Header: string "SWEETFileDict" which is 0 terminated
    
    Number of entries: 64 bit integer
    
    For number of entries:
    
        Key string: string which is 0 terminated
        
        type_id:
            100: string which is 0 terminated
            230: scalar 64 bit integer (64 bit)
            240: scalar 64 bit double float (64 bit)
            401: 1D array with floating point numbers (each 64 bit)
            402: 2D array with floating point numbers (each 64 bit)
            403: 3D array with floating point numbers (each 64 bit)
            501: 1D array with complex numbers (2x 64 bit)
            502: 2D array with complex numbers (2x 64 bit)
            503: 3D array with complex numbers (2x 64 bit)
            
        For type_id >= 200:
            This is followed by a number of 64 bit integers which describe the size in each dimension
            
        After this, the raw data is written.
    
    Final magic code: string "SWEETFileDict" with 0 terminal
    """
    def __init__(self, filename=None, debug=False, initDict=None):
        self.dict = {}
        if initDict is not None:
            self.dict.update(initDict)
        
        self.debug = debug
        
        if filename != None:
            self.readFromFile(filename)

    def update(self, dict):
        """Update the SWEETFileDict with a given Python dictionnary"""
        self.dict.update(dict)
        
    def set(self, key, value):
        """Update (or add) one entry of the SWEETFileDict"""
        self.dict[key] = value
        
    def delete(self, key):
        del self.dict[key]

    def __setitem__(self, key, value):
        self.set(key, value)

    def __getitem__(self, key):
        return self.dict[key]
    
    def __str__(self):
        retstr = ""
        for key in self.dict:
            retstr += f"'{key}' => '{self.dict[key]}'\n"

        return retstr

    def writeToFile(self, filename):
        f = open(filename, "wb")
        
        self._write_str0(f, "SWEETFileDict")
        
        self._write_int64(f, len(self.dict));
        
        for key in self.dict:
            
            # Write key
            self._write_str0(f, key)
            
            value = self.dict[key]
            
            if isinstance(value, str):
                self._write_int64(f, 100)
                self._write_str0(f, value)
                
            elif isinstance(value, int):
                self._write_int64(f, 230)
                self._write_int64(f, value)
            
            elif isinstance(value, float):
                self._write_int64(f, 240)
                self._write_double(f, value)

            elif isinstance(value, np.ndarray):
                if value.dtype == np.float64:
                    if value.ndim == 1:
                        self._write_int64(f, 401)
                    elif value.ndim == 2:
                        self._write_int64(f, 402)
                    elif value.ndim == 3:
                        self._write_int64(f, 403)
                    else:
                        raise Exception(f"Key '{key}': Unsupported type {type(value)} of dimension {value.ndim}")

                    for i in range(value.ndim):
                        self._write_int64(f, value.shape[i])
                    
                    self._write_ndarray_double(f, value)
                    
                elif value.dtype == np.complex128:
                    if value.ndim == 1:
                        self._write_int64(f, 501)
                    elif value.ndim == 2:
                        self._write_int64(f, 502)
                    elif value.ndim == 3:
                        self._write_int64(f, 503)
                    else:
                        raise Exception(f"Key '{key}': Unsupported type {type(value)} of dimension {value.ndim}")

                    for i in range(value.ndim):
                        self._write_int64(f, value.shape[i])
                    
                    self._write_ndarray_complex(f, value)
                    
                else:
                    raise Exception(f"Key '{key}': Unsupported type {np.dtype} of ndarray")
                

            else:
                raise Exception(f"Key '{key}': Unsupported type {type(value)}")
        
            if self.debug:
                print(f"Written {key}")
            
        self._write_str0(f, "SWEETFileDict")

    def _write_byte(self, f, value):
        f.write(struct.pack('b', value))
                
    def _write_str(self, f, string):
        f.write(string.encode('ascii'))
        
    def _write_str0(self, f, string):
        self._write_str(f, string)
        self._write_byte(f, 0)
        
    def _write_int64(self, f, value):
        f.write(struct.pack('q', value))
        
    def _write_double(self, f, value):
        f.write(struct.pack('d', value))

    def _write_ndarray_double(self, f, value):
        v = value.flatten()
        for i in range(v.size):
            self._write_double(f, v[i])

    def _write_complex(self, f, value):
        self._write_double(f, value.real)
        self._write_double(f, value.imag)
        

    def _write_ndarray_complex(self, f, value):
        v = value.flatten()
        for i in range(v.size):
            self._write_complex(f, v[i])
            
            
    def readFromFile(self, filename):
        f = open(filename, "rb")
        
        # Read Magic code header
        value = self._read_str(f)
        assert value == "SWEETFileDict", "Magic code not found!"
        
        num_keys = self._read_int64(f)
        
        if self.debug:
            print(f" + keys: {num_keys}")
        
        for key_id in range(num_keys):
            
            # read key
            key = self._read_str(f)

            type_id = self._read_int64(f)    

            if self.debug:
                print(f" + key: {key} with type_id {type_id}")
        
            if type_id == 100:
                # string
                self.dict[key] = self._read_str(f)
                
            elif type_id == 230:
                self.dict[key] = self._read_int64(f)

            elif type_id == 240:
                self.dict[key] = self._read_double(f)
                
            elif type_id == 401:
                # 1D np.array of type double
                shape = (self._read_int64(f), )
                self.dict[key] = self._read_ndarray_double(f, shape)
                
            elif type_id == 402:
                # 2D np.array of type double
                shape = (self._read_int64(f), self._read_int64(f), )
                self.dict[key] = self._read_ndarray_double(f, shape)
                
            elif type_id == 403:
                # 3D np.array of type double
                shape = (self._read_int64(f), self._read_int64(f), self._read_int64(f), )
                self.dict[key] = self._read_ndarray_double(f, shape)
                
                
            elif type_id == 501:
                # 1D np.array of type complex
                shape = (self._read_int64(f), )
                self.dict[key] = self._read_ndarray_complex(f, shape)
                
            elif type_id == 502:
                # 2D np.array of type complex
                shape = (self._read_int64(f), self._read_int64(f), )
                self.dict[key] = self._read_ndarray_complex(f, shape)
                
            elif type_id == 503:
                # 3D np.array of type complex
                shape = (self._read_int64(f), self._read_int64(f), self._read_int64(f), )
                self.dict[key] = self._read_ndarray_complex(f, shape)

            else:
                raise Exception(f"type id '{type_id}' not supported")

            if self.debug:
                print(f" ++ Read {key} with value '{self.dict[key]}'")
            
        fin_str = self._read_str(f)
        
        assert fin_str == "SWEETFileDict", "Final magic code not found!"

    def _read_byte(self, f, value):
        return struct.unpack('b', f.read(1))[0]
                
    def _read_str(self, f):
        s = bytearray()
        while True:
            b = f.read(1)
            if b == b'\x00':
                break
            s.extend(b)
        return s.decode()

    def _read_int64(self, f):
        return struct.unpack('q', f.read(8))[0]
        
    def _read_double(self, f):
        return struct.unpack('d', f.read(8))[0]

    def _read_ndarray_double(self, f, shape):    
        size = np.prod(shape)
        fa = np.empty(dtype=np.float64, shape=size)
        
        for i in range(size):
            fa[i] = self._read_double(f)

        return fa.reshape(shape)
    
    def _read_complex(self, f):
        return struct.unpack('d', f.read(8))[0]+struct.unpack('d', f.read(8))[0]*1j
        

    def _read_ndarray_complex(self, f, shape):
        size = np.prod(shape)
        fa = np.empty(dtype=np.complex128, shape=size)
        
        for i in range(size):
            fa[i] = self._read_complex(f)

        return fa.reshape(shape)

    def __eq__(self, other):
        if len(self.dict) != len(other.dict):
            return False
        
        for key in self.dict:
            
            value = self.dict[key]
            
            if not isinstance(value, np.ndarray):
                if self.dict[key] != other.dict[key]:
                    return False
            else:
                if np.any(self.dict[key] != other.dict[key]):
                    return False
        
        return True