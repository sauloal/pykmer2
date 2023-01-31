#!/usr/bin/env python3

#https://github.com/sauloal/pykmer2

import os
import sys
import math
import datetime
import struct
import functools

"""
1
4.0K idx_03
4.0K idx_03.xz

2
4.0K idx_05
4.0K idx_05.xz

2
 12K idx_07
8.0K idx_07.xz

4
260K idx_09
 28K idx_09.xz

4
4.1M idx_11
284K idx_11.xz

4
 65M idx_13
4.1M idx_13.xz

4
1.1G idx_15
 61M idx_15.xz

8
idx_17
"""

"""
time pypy3 kmer-numbering.py 11
 4**11    =      4,194,304
 4**11/4  =      1,048,576
 4**11/4*4=      4,194,304b 4,09Kkb 4Mb
 num_regs =      1,049,600
 i        =      4,100,000
real    2m19.662s
user    2m18.588s
sys     0m 1.031s

time pypy3 kmer-numbering.py 15
 4**15    =   1,073,741,824
 4**15/4  =     268,435,456
 4**15/4*4=   1,073,741,824b 1,048,576Kb 1,024Mb 1Gb
 num_regs =     268,451,840
 i        =   1,073,700,000
real    1703m 3.469s
user    1575m32.661s
sys       66m47.063s

time pypy3 kmer-numbering.py 17
 4**17    = 17,179,869,184
 4**17/4  =  4,294,967,296
 4**17/4*8= 34,359,738,368b 33,554,432Kb 32,768Mb 32Gb
 num_regs = 

time pypy3 kmer-numbering.py 19
 4**19    = 274,877,906,944
 4**19/4  =  68,719,476,736
 4**19/4*8= 549,755,813,888b 536870912Kb 524,288Mb 512Gb
 num_regs = 
 i        = 

time pypy3 kmer-numbering.py 21
 4**21    = 
 4**21/4  = 
 4**21/4*8= 8,589,934,592Kb 8,388,608Mb 8,192Gb 8Tb
 num_regs = 
 i        = 


"""

VOCAB = ['A','C','G','T']
RCD   = {'A':'T','C':'G','G':'C','T':'A'}
RC    = [RCD.get(chr(c), None) for c in range(254)]


BACOVD = {     c:i for i,c in enumerate(VOCAB) }
BACOV  = [BACOVD.get(chr(c), None) for c in range(254)]

MULTS = tuple(i   for i   in range(32,-2,-2))
STLUM = tuple(reversed(MULTS))
#print("MULTS", MULTS)

MATRIX = [None] * 255
for i,c in enumerate(VOCAB):
	ordc = ord(c)
	print(f"{i=} {c=} {ordc=}")
	MATRIX[ordc] = [None] * 32
	for pos in range(32):
		v = (4**(32-pos-1)) * i
		MATRIX[ordc][pos] = v
#print("MATRIX", MATRIX)

struct_formats = {
	 3: ">B", # log2(4** 3) =  6 1
	 5: ">H", # log2(4** 5) = 10 2
	 7: ">H", # log2(4** 7) = 14 2
	 9: ">L", # log2(4** 9) = 18 4
	11: ">L", # log2(4**11) = 22 4
	13: ">L", # log2(4**13) = 26 4
	15: ">L", # log2(4**15) = 30 4
	17: ">Q", # log2(4**17) = 34 8
	19: ">Q", # log2(4**19) = 38 8
	21: ">Q", # log2(4**21) = 42 8
	23: ">Q", # log2(4**23) = 46 8
	25: ">Q", # log2(4**25) = 50 8
	27: ">Q", # log2(4**27) = 54 8
	29: ">Q", # log2(4**29) = 58 8
	31: ">Q", # log2(4**31) = 62 8
}

def generate_kmer(length, pos=0, prev=None):
	if pos == length:
		#print("length", prev)
		yield prev
		return

	for c in VOCAB:
		#print("pos", pos, "c", c)
		for p in generate_kmer(length, pos=pos+1, prev=("" if prev is None else prev)+c):
			#print(" p", p)
			yield p

#@functools.lru_cache(maxsize=None)
def get_mask(kmer_size):
	masks = [3 << ((kmer_size-k-1)*2) for k in range(kmer_size)]
	return masks

def generate_sequence(kmer_size, value):
	#TODO: cache each byte
	masks = get_mask(kmer_size)
	#print(masks, [f"{m:>06b}" for m in masks])
	nucs = [None]*kmer_size
	for k in range(kmer_size):
		mask    = masks[k]
		val     = value & mask
		nic     = val >> ((kmer_size-k-1)*2)
		nuc     = VOCAB[nic]
		nucs[k] = nuc
		#print(f"{value=:02d} {value:>06b} {k=} {mask=:02d} {mask:>06b} {val=:02d} {val:>06b} {nic=:02d} {nic:>06b} {nuc} {nucs}")
	return "".join(nucs)

@functools.lru_cache(maxsize=1_000_000)
def index_kmer(seq):
	lseq     = len(seq)
	mults    = STLUM[:lseq]
	muls     = tuple(mults[lseq-i-1] for i in range(lseq))

	values   = (BACOV[c]    for c   in seq)
	calcs    = (v<<muls[i]  for i,v in enumerate(values))
	calc_sum = sum(calcs)
	#print("  ", seq, values, calcs, calc_sum)
	return calc_sum

def get_kmer_indexer(lseq):
	mults    = STLUM[:lseq]
	muls     = tuple(mults[lseq-i-1] for i in range(lseq))

	#@functools.lru_cache(maxsize=1_000_000)
	def kmer_indexer(seq):
		values   = (BACOV[ord(c)] for c   in seq)
		calcs    = (v<<muls[i]    for i,v in enumerate(values))
		calc_sum = sum(calcs)
		#print("  ", seq, values, calcs, calc_sum)
		return calc_sum

	return kmer_indexer

kmer_indexer = None

#@functools.lru_cache(maxsize=1_000_000)
def rev_comp_4(seq, debug=False):
	fwd      = tuple(seq)
	fwd_comp = tuple(RC[ord(c)] for c in fwd)
	rev      = tuple(fwd[::-1])
	rev_comp = tuple(fwd_comp[::-1])

	#mask         = (2 << ((len(seq)-1)*2))
	#print(f"  MASK    {mask:5d} {mask:>06b} {mask:03x}")

	global kmer_indexer
	if kmer_indexer is None:
		kmer_indexer = get_kmer_indexer(len(fwd))

	fwd_idx      = kmer_indexer(fwd)
	fwd_comp_idx = kmer_indexer(fwd_comp)
	rev_idx      = kmer_indexer(rev)
	rev_comp_idx = kmer_indexer(rev_comp)

	is_fwd  = None
	is_comp = None
	is_fake = None
	min_seq = None
	min_idx = None

	if   fwd_idx      <= fwd_comp_idx and fwd_idx      <= rev_idx      and fwd_idx      <= rev_comp_idx: # fwd
		is_fwd  = True
		is_comp = False
		is_fake = False
		min_seq = fwd
		min_idx = fwd_idx
	elif rev_comp_idx <=  fwd_idx     and rev_comp_idx <= fwd_comp_idx and rev_comp_idx <= rev_idx: # rev_comp
		is_fwd  = False
		is_comp = True
		is_fake = False
		min_seq = rev_comp
		min_idx = rev_comp_idx
	elif fwd_comp_idx <= fwd_idx      and fwd_comp_idx <= rev_idx      and fwd_comp_idx <= rev_comp_idx: # fwd_comp
		is_fwd  = True
		is_comp = True
		is_fake = True
		min_seq = fwd_comp
		min_idx = fwd_comp_idx
	elif rev_idx      <= fwd_idx      and rev_idx      <= fwd_comp_idx and rev_idx      <= rev_comp_idx: # rev
		is_fwd  = False
		is_comp = False
		is_fake = True
		min_seq = rev
		min_idx = rev_idx

	min_seq = "".join(min_seq)

	if debug: print(f"  FWD      {fwd     =} {fwd_idx     =:5d} {fwd_idx     :>06b} {fwd_idx     :03x} {'*' if     is_fwd and not is_comp else ''}")
	if debug: print(f"  FWD COMP {fwd_comp=} {fwd_comp_idx=:5d} {fwd_comp_idx:>06b} {fwd_comp_idx:03x} {'*' if     is_fwd and     is_comp else ''}")
	if debug: print(f"  REV      {rev     =} {rev_idx     =:5d} {rev_idx     :>06b} {rev_idx     :03x} {'*' if not is_fwd and not is_comp else ''}")
	if debug: print(f"  REV COMP {rev_comp=} {rev_comp_idx=:5d} {rev_comp_idx:>06b} {rev_comp_idx:03x} {'*' if not is_fwd and     is_comp else ''}")

	return min_seq, min_idx, is_fwd, is_comp, is_fake

def calc_index_algo(kmer_size, seq, debug=False):
	debug     = True
	fwd       = tuple(seq)
	rec       = tuple(RC[c] for c in fwd[::-1])
	fwd_idx   = index_kmer(fwd)
	rec_idx   = index_kmer(rec)

	ceq,idx   = (fwd,fwd_idx) if fwd_idx <= rec_idx else (rec,rec_idx)
	idx_block = idx // 4
	idx_mod   = idx %  4

	if debug: print(f"    {kmer_size=} {fwd=} {fwd_idx=} {rec=} {rec_idx=} {ceq=} {idx=}")

	return calc_index_algo_idx(kmer_size, idx, debug=debug)

def calc_index_algo_idx(kmer_size, idx, debug=False):
	debug               = True

	idx_col             = idx // 4
	idx_mod             = idx %  4

	if debug: print(f"    {idx=} {idx_col=} {idx_mod=}")

	mod_block_cols        = 4
	mod_block_rows        = 4
	mod_block_size        = mod_block_cols * mod_block_rows

	mod_block_offset      = 2
	mod_block_offset_cols = mod_block_offset * mod_block_cols

	mod_offset_1          = idx          - mod_block_offset_cols + 1
	mod_offset_2          = mod_offset_1 - mod_block_size        + 1
	mod_offset_3          = mod_offset_2 - mod_block_size        + 1
	mod_offset_4          = mod_offset_3 - mod_block_size        + 1

	mod_offset_1_mod      = mod_offset_1 // 4
	mod_offset_2_mod      = mod_offset_2 // 4
	mod_offset_3_mod      = mod_offset_3 // 4
	mod_offset_4_mod      = mod_offset_4 // 4

	if debug: print(f"      {mod_block_cols=:2d} {mod_block_rows=:2d} {mod_block_size=:2d} {mod_block_offset=:2d}")
	if debug: print(f"        {mod_offset_1=:5d} {mod_offset_1_mod=:5d}")
	if debug: print(f"        {mod_offset_2=:5d} {mod_offset_2_mod=:5d}")
	if debug: print(f"        {mod_offset_3=:5d} {mod_offset_3_mod=:5d}")
	if debug: print(f"        {mod_offset_4=:5d} {mod_offset_4_mod=:5d}")

	idx_delta = max(0,mod_offset_1_mod) + max(0,mod_offset_2_mod) + max(0,mod_offset_3_mod) + max(0,mod_offset_4_mod)
	idx_end   = idx - max(0,idx_delta)
	if debug: print(f"    {idx_delta=:2d} {idx_end=:2d}")
	if debug: print()

	"""
	boundary_mod_offset_global = 2
	boundary_mod_offset_block  = 3
	boundary_mod_block_size    = 4
	boundary_mod_pilar_size    = 4**(kmer_size-3)
	boundaries                 = gen_boundaries(kmer_size)

	#if debug: print(f"  {boundary_mod_offset_global=} {boundary_mod_offset_block=} {boundary_mod_block_size=} {boundary_mod_pilar_size=}")
	#if debug: print(f"  {boundaries=}")

	idx_mods = [0] * boundary_mod_block_size
	for i in range(boundary_mod_block_size):
		#=(idx_mod+1)*if(idx_mod>boundary_pos,1,-1)
		if False:
			boundary_pos         = boundary_mod_offset_global + (boundary_mod_block_size*i)
			boundary_col         = boundary_mod_block_size - i
			boundary_pos_valid   =  idx_block  >= boundary_pos
			boundary_col_valid   = (idx_mod+1) == boundary_col
			boundary_block = -1 # ((boundary_mod_offset_global + (i*4))*4) - i
			boundary       = -1 # boundary_mod_offset_block  + boundary_block
			diff_mod       = (idx - boundary) // 4 if idx_mod >= boundary else 0
			idx_mods[i]    = diff_mod
			if debug: print(f"  {i=:2d} {boundary_pos=:2d} {boundary_col=:2d} {boundary_pos_valid=:2d} {boundary_col_valid=:2d} {boundary_block=:2d} {boundary=:2d} {diff_mod=:2d}")

	idx_sum = sum(idx_mods)
	idx_end = idx - idx_sum

	print(f"  {idx_sum=} {idx_end=}")
	"""

	return idx_end

def gen_boundaries(kmer_size, size=4, debug=False):
	block_size                   = 4
	boundaries_mod_offset_global = 2
	boundaries_mod_offset_block  = 3
	boundaries_mod               = [ boundaries_mod_offset_block+((boundaries_mod_offset_global+(v*4))*4)-v for v in range(size) ]
	if debug: print("  boundaries_mod", boundaries_mod)
	return boundaries_mod

def calc_index(ind, kmer_size, debug=False):
	boundaries_mod = gen_boundaries(kmer_size, debug=debug)

	if debug: print("\n ind", ind)
	diff_mod                 = [ (ind-b+4)//4 if ind >= b else 0 for i,b in enumerate(boundaries_mod) ]
	if debug: print("  diff_mod      ", diff_mod)
	idx  = ind - sum(diff_mod)
	if debug: print("  idx", idx)
	return idx

def calc_index_2(seq, idx, kmer_size, debug=False):
	#debug = True

	if debug: print()
	if debug: print("calc_index_2", seq)

	sen = [BACOV[s] for s in seq]
	val = [0 for s in sen]

	if debug: print(f"  {seq=} {idx=}")
	if debug: print(f"  idx  {idx :12d} {idx :>06b}")

	if idx < 32:
		#first letter is always a or c
		mask = (1 << ((kmer_size) * 2)-1)-1 # clear first bit 011111
		val  = idx & mask
		if debug: print(f"  mask {mask:12d} {mask:>06b}")
		if debug: print(f"  val  {val :12d} {val :>06b}")
	else:
		#last char is always a or c

		mask =  idx & 1 # get last bit. 0 for A, 1 for C - 000001
		lval = (idx & mask)
		val  = (idx >> 1) | lval # delete last bit, move to last-to-last bit, converting T to A and G to C
		if debug: print(f"  mask {mask:12d} {mask:>06b}")
		if debug: print(f"  lval {lval:12d} {lval:>06b}")
		if debug: print(f"  val  {val :12d} {val :>06b}")

	return val

#def rev_comp(seq):
#	return "".join([RC[c] for c in seq[::-1]])

class Kmer:
	def __init__(self, kmer_size, debug=False):
		self.kmer_size     = kmer_size
		self.debug         = debug

		self.struct_fmt    = struct_formats[self.kmer_size]
		self.struct_size   = struct.calcsize(self.struct_fmt)
		self.struct_pack   = struct.Struct(self.struct_fmt).pack
		self.struct_unpack = struct.Struct(self.struct_fmt).unpack_from
		self.struct_max    = 2**(self.struct_size*8)-1

		print(f"{self.struct_fmt=} {self.struct_size=} {self.struct_max=:15,d}")

		self.index_file = f"idx_{kmer_size:02d}"
		self.fhd_w = None
		self.fhd_r = None
		self.num_regs = 0

	@property
	def data_pos(self):
		return (3*8)
	@property
	def footer_pos(self):
		return (3*8) + (self.num_regs*self.struct_size)

	@property
	def exists(self):
		return os.path.exists(self.index_file)

	def open_w(self):
		assert self.fhd_w is None
		print("opening")
		self.fhd_w = open(self.index_file+'.tmp', "wb")

		self.fhd_w.write(struct.pack('>Q',self.struct_max))
		self.fhd_w.write(struct.pack('>Q',self.kmer_size))
		self.fhd_w.write(struct.pack('>Q',0))

		self.num_regs = 0

	def close_w(self):
		assert self.fhd_w is not None
		print("closing")
		self.fhd_w.write(struct.pack('>Q',self.num_regs))
		self.fhd_w.write(struct.pack('>Q',self.kmer_size))
		self.fhd_w.write(struct.pack('>Q',self.struct_max))
		self.fhd_w.flush()
		self.fhd_w.seek(2*8)
		self.fhd_w.write(struct.pack('>Q',self.num_regs))
		self.fhd_w.flush()
		self.fhd_w.close()
		self.fhd_w = None

		os.rename(self.index_file+'.tmp',self.index_file)
		print("closed")

	def pack_w(self, val):
		assert self.fhd_w is not None
		assert val <= self.struct_max
		self.num_regs += 1
		self.fhd_w.write(self.struct_pack(val))

	def open_r(self, list_regs=False, kmer_size=None, num_regs=None):
		assert self.fhd_w is None
		assert self.fhd_r is None

		self.fhd_r       = open(self.index_file, "rb")

		filesize         = os.path.getsize(self.index_file)

		struct_max_val_h = struct.unpack('>Q', self.fhd_r.read(8))[0]
		kmer_size_val_h  = struct.unpack('>Q', self.fhd_r.read(8))[0]
		num_regs_val_h   = struct.unpack('>Q', self.fhd_r.read(8))[0]

		assert struct_max_val_h == self.struct_max
		assert kmer_size  is None or kmer_size_val_h  == kmer_size , f"{kmer_size_val_h=} {kmer_size=}"
		assert num_regs   is None or num_regs_val_h   == num_regs  , f"{num_regs_val_h=} {num_regs=}"

		print(f"{struct_max_val_h=:15,d}")
		print(f"{kmer_size_val_h=:15,d}")
		print(f"{num_regs_val_h=:15,d}")

		num_regs      = (filesize - (6*8)) // self.struct_size
		assert num_regs == num_regs_val_h
		self.num_regs = num_regs

		if list_regs:
			self.list_regs()
			assert self.fhd_r.tell() == self.footer_pos

		self.fhd_r.seek(self.footer_pos)

		num_regs_val_t   = struct.unpack('>Q', self.fhd_r.read(8))[0]
		kmer_size_val_t  = struct.unpack('>Q', self.fhd_r.read(8))[0]
		struct_max_val_t = struct.unpack('>Q', self.fhd_r.read(8))[0]

		assert num_regs   is None or num_regs_val_t   == num_regs  , f"{num_regs_val_t=} {num_regs=}"
		assert kmer_size  is None or kmer_size_val_t  == kmer_size , f"{kmer_size_val_t=} {kmer_size=}"

		assert struct_max_val_t == self.struct_max
		assert num_regs_val_t   == num_regs_val_h
		assert kmer_size_val_t  == kmer_size_val_h
		assert struct_max_val_t == struct_max_val_h

		print("valid")

		self.fhd_r.seek(self.data_pos)

	def close_r(self):
		assert self.fhd_r is not None
		print("closing")
		self.fhd_r.close()
		self.fhd_r = None

	def list_regs(self):
		assert self.fhd_r is not None

		for reg_num in range(self.num_regs):
			reg = self.fhd_r.read(self.struct_size)
			if not reg: break
			val = self.struct_unpack(reg)[0]
			if self.debug or reg_num % 1_000 == 0:
				print(f"{reg_num=:15,d} {str(reg)=:4s} {val=:15,d}")

		assert self.fhd_r.tell() == self.footer_pos


	#@functools.lru_cache(maxsize=1_000_000)
	def get_register_at_pos(self, pos):
		assert self.fhd_r is not None
		assert pos < self.num_regs, f"{pos=} < {self.num_regs=}"
		relapos = pos*self.struct_size
		self.fhd_r.seek(self.data_pos + relapos)
		reg = self.fhd_r.read(self.struct_size)
		val = self.struct_unpack(reg)[0]
		return val

	#@functools.lru_cache(maxsize=1_000_000)
	def find_register_id(self, value):
		lo = 0
		hi = self.num_regs
		while True:
			#print(f"{lo=} {hi=} {value=}")
			if lo == hi:
				return None, None
			mi = (hi + lo) // 2
			re = self.get_register_at_pos(mi)
			#print(f"  {mi=} {re=}")
			if   re == value:
				return mi, re
			elif re > value:
				hi = mi
			else:
				lo = mi

	#@functools.lru_cache(maxsize=1_000_000)
	def find_sequence_id(self, seq):
		assert self.fhd_r is not None
		cds, cdx, is_fwd, is_comp, is_fake = rev_comp_4(seq, debug=False)
		pos, val                           = self.find_register_id(cdx)
		qes                                = self.generate_sequence(cdx)
		assert val == cdx
		assert qes == cds
		if self.debug: print(f"  {cds} {cdx:{int(math.log10(self.struct_max))}d} {cdx:>0{self.struct_size*8}b} {is_fwd=!s:6s} {is_comp=!s:6s} {is_fake=!s:6s} {pos=:{int(math.log10(self.struct_max))}d} {val=:{int(math.log10(self.struct_max))}d}")
		return pos, val

	#@functools.lru_cache(maxsize=1_000_000)
	def generate_sequence(self, value):
		return generate_sequence(self.kmer_size, value)

	#@functools.lru_cache(maxsize=1_000_000)
	def rev_comp_4(self, seq):
		return rev_comp_4(seq, debug=self.debug)

def create2(kmer, debug=False):
	for i, seq in enumerate(generate_kmer(kmer.kmer_size)):
		calc_index_algo_idx(kmer_size, i, debug=debug)
		#calc_index_algo(kmer.kmer_size, seq, debug=debug)

class Eta:
	def __init__(self, total=-1):
		self._start  = datetime.datetime.now()
		self._last_t = self._start
		self._last_v = -1
		self._total  = total

		self.now           = datetime.datetime.now()
		self.delta_t_start = datetime.timedelta(seconds=-1)
		self.delta_t_last  = datetime.timedelta(seconds=-1)
		self.diff_v_start  = -1
		self.diff_v_last   = -1
		self.speed_start   = -1
		self.speed_last    = -1
		self.eta_v         = -1
		self.eta_t         = -1

	def update(self, value):
		self.now           = datetime.datetime.now()
		self.delta_t_start = self.now - self._start
		self.delta_t_last  = self.now - self._last_t

		if value < self._last_v:
			self.diff_v_last = -1
			self.speed_last  = -1
		else:
			ela_last         = self.delta_t_last.total_seconds()
			self.diff_v_last = value - self._last_v
			self.speed_last  = self.diff_v_last / ela_last


		if   self._total == -1:
			self.eta_v = -1
			self.eta_t = -1
		else:
			ela_start          = self.delta_t_start.total_seconds()
			self.diff_v_start  = value
			self.speed_start   = self.diff_v_start / ela_start
			self.eta_v         = self._total - value

			if self.speed_start == 0:
				self.eta_t = -1
			else:
				etas               = self.eta_v / self.speed_start
				self.eta_t         = datetime.timedelta(seconds=int(etas))

		self._last_t = self.now
		self._last_v = value

	def __str__(self):
		now           = self.now.isoformat(sep=" ", timespec="seconds")
		delta_t_start = datetime.timedelta(days=self.delta_t_start.days, seconds=self.delta_t_start.seconds)
		delta_t_last  = datetime.timedelta(days=self.delta_t_last.days , seconds=self.delta_t_last.seconds)
		return f"{now} | {self._last_v:15,d} | Delta V :: Start/Last {self.diff_v_start:15,d}/{self.diff_v_last:15,d} | Delta T: Start/Last {str(delta_t_start):15s}/{str(delta_t_last):15s} | Speed Start/Last {int(self.speed_start):15,d}/{int(self.speed_last):15,d} | ETA {self.eta_t}/{self.eta_v:15,d}"

def create(kmer, debug=False):
	print("creating")

	ids = {}
	kmer.open_w()

	eta = Eta(4**kmer_size)

	num_regs = 0
	for i, seq in enumerate(generate_kmer(kmer.kmer_size)):
		if False:
			qes = rev_comp(seq)
			ind = index_kmer(seq)
			qes = rev_comp(seq)
			dni = index_kmer(qes)
			inv = ind > dni
			sex = seq if not inv else qes
			idx = ind if not inv else dni
			#cdx = calc_index(idx, kmer_size, debug=False)
			cdx = calc_index_2(sex, idx, kmer_size)
			if not inv:
				if debug: print(f"SEQ {i=:5d} {seq=} {ind=:5d} {ind:>06b} {inv=!s:6s} {qes=} {dni=:5d} {dni:>06b} {qes if inv else seq} {idx=:5d} {idx:>06b} {idx:03x} {cdx=:5d} {cdx:>06b} {cdx:03x} {'*' if inv else ''}")
		else:
			if debug: print(f"SEQ {i=:5d} {seq=}")
			cds, cdx, is_fwd, is_comp, is_fake = kmer.rev_comp_4(seq)
			if debug: print(f"  {cds} {cdx:5d} {cdx:>06b} {cdx:03x} {is_fwd=!s:6s} {is_comp=!s:6s} {is_fake=!s:6s}")

		if kmer_size<7:
			ids[cdx] = ids.get(cdx, []) + [i]

		if i % (100_000 if (4**kmer_size) > 100_000 else 1) == 0:
			print(f" {eta} | {num_regs=:15,d}")
			eta.update(i)

		if cdx >= i:
			num_regs += 1
			kmer.pack_w(cdx)

	print(f"{num_regs=:15,d}")

	kmer.close_w()

	print("created")

	if kmer.kmer_size < 7:
		for i, (cdx, idxs) in enumerate(sorted(ids.items())):
			print(f"{i:3d} {cdx:3d} {idxs}")

	return num_regs, ids

def check(kmer, num_regs=None, debug=False):
	kmer.open_r(list_regs=kmer.kmer_size<7, kmer_size=kmer.kmer_size, num_regs=num_regs)

	#kmer.list(kmer_size=kmer_size, list_regs=kmer_size<7, debug=kmer_size<7)

	print("checking")
	for i, seq in enumerate(generate_kmer(kmer.kmer_size)):
		if i % (100_000 if (4**kmer_size) > 100_000 else 1) == 0:
			print(f" {i=:15,d} / {4**kmer_size:15,d}")
		pos, val = kmer.find_sequence_id(seq)
		#calc_index_algo(kmer.kmer_size, seq, debug=debug)
		#cds, cdx, is_fwd, is_comp, is_fake = kmer.rev_comp_4(seq)
		#qes                                = kmer.generate_sequence(cdx)
		#if debug or kmer_size<7: print(f"{i=} {seq=} {pos=} {val=} {cds=} {cdx=} {qes=}")
		#assert qes == cds, f'{qes=} == {cds=}'
	print("checked")

	kmer.close_r()

def main(kmer_size, debug=False):
	kmer = Kmer(kmer_size=kmer_size, debug=kmer_size<7 or debug)

	num_regs = None


	#create2(kmer, debug=debug)
	#sys.exit(0)

	if kmer.exists:
		print("exists")
	else:
		num_regs, ids = create(kmer, debug=debug)


	if kmer_size < 15:
		#check(kmer, num_regs=num_regs, debug=debug)
		pass
	else:
		print("too large. not checking")

if __name__ == "__main__":
	kmer_size = int(sys.argv[1])
	main(kmer_size=kmer_size)

