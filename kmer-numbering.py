#!/usr/bin/env python3

import os
import sys
import math
import struct
import functools

"""
4.0K idx_03
4.0K idx_03.xz

4.0K idx_05
4.0K idx_05.xz

 12K idx_07
8.0K idx_07.xz

260K idx_09
 28K idx_09.xz

4.1M idx_11
284K idx_11.xz

 65M idx_13
4.1M idx_13.xz

1.1G idx_15
 61M idx_15.xz
"""

"""
time pypy3 kmer-numbering.py 11
 num_regs=   1,049,600
 i       =   4,100,000
real    2m19.662s
user    2m18.588s
sys     0m 1.031s

time pypy3 kmer-numbering.py 15
 num_regs= 268,451,840

"""

VOCAB = ['A','C','G','T']
RC    = {'A':'T','C':'G','G':'C','T':'A'}

BACOV = { c:i for i,c in enumerate(VOCAB) }
MULTS = tuple(i      for i in range(32,-2,-2))
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

@functools.lru_cache(maxsize=None)
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
	values   = [BACOV[c] for c in seq]
	calcs    = [v<<mults[lseq-i-1] for i,v in enumerate(values)]
	calc_sum = sum(calcs)
	#print("  ", seq, values, calcs, calc_sum)
	return calc_sum

@functools.lru_cache(maxsize=1_000_000)
def rev_comp_4(seq, debug=False):
	fwd      = tuple(seq)
	fwd_comp = tuple(RC[c] for c in fwd)
	rev      = tuple(fwd[::-1])
	rev_comp = tuple(fwd_comp[::-1])

	#mask         = (2 << ((len(seq)-1)*2))
	#print(f"  MASK    {mask:5d} {mask:>06b} {mask:03x}")

	fwd_idx      = index_kmer(fwd)
	fwd_comp_idx = index_kmer(fwd_comp)
	rev_idx      = index_kmer(rev)
	rev_comp_idx = index_kmer(rev_comp)

	is_fwd  = None
	is_comp = None
	is_fake = None
	min_seq = None
	min_idx = None

	if   all(fwd_idx      <= r for r in (fwd_comp_idx, rev_idx     , rev_comp_idx)): # fwd
		is_fwd  = True
		is_comp = False
		is_fake = False
		min_seq = fwd
		min_idx = fwd_idx
	elif all(rev_comp_idx <= r for r in (fwd_idx     , fwd_comp_idx, rev_idx     )): # rev_comp
		is_fwd  = False
		is_comp = True
		is_fake = False
		min_seq = rev_comp
		min_idx = rev_comp_idx
	elif all(fwd_comp_idx <= r for r in (fwd_idx     , rev_idx     , rev_comp_idx)): # fwd_comp
		is_fwd  = True
		is_comp = True
		is_fake = True
		min_seq = fwd_comp
		min_idx = fwd_comp_idx
	elif all(rev_idx      <= r for r in (fwd_idx     , fwd_comp_idx, rev_comp_idx)): # rev
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


	@functools.lru_cache(maxsize=1_000_000)
	def get_register_at_pos(self, pos):
		assert self.fhd_r is not None
		assert pos < self.num_regs, f"{pos=} < {self.num_regs=}"
		relapos = pos*self.struct_size
		self.fhd_r.seek(self.data_pos + relapos)
		reg = self.fhd_r.read(self.struct_size)
		val = self.struct_unpack(reg)[0]
		return val

	@functools.lru_cache(maxsize=1_000_000)
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

	@functools.lru_cache(maxsize=1_000_000)
	def find_sequence_id(self, seq):
		assert self.fhd_r is not None
		cds, cdx, is_fwd, is_comp, is_fake = rev_comp_4(seq, debug=False)
		pos, val                           = self.find_register_id(cdx)
		qes                                = self.generate_sequence(cdx)
		assert val == cdx
		assert qes == cds
		if self.debug: print(f"  {cds} {cdx:{int(math.log10(self.struct_max))}d} {cdx:>0{self.struct_size*8}b} {is_fwd=!s:6s} {is_comp=!s:6s} {is_fake=!s:6s} {pos=:{int(math.log10(self.struct_max))}d} {val=:{int(math.log10(self.struct_max))}d}")
		return pos, val

	@functools.lru_cache(maxsize=1_000_000)
	def generate_sequence(self, value):
		return generate_sequence(self.kmer_size, value)

	@functools.lru_cache(maxsize=1_000_000)
	def rev_comp_4(self, seq):
		return rev_comp_4(seq, debug=self.debug)


def create(kmer, debug=False):
	print("creating")

	ids = {}
	kmer.open_w()

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
			cds, cdx, is_fwd, is_comp, is_fake = kmer.rev_comp_4(seq, debug=False)
			if debug: print(f"  {cds} {cdx:5d} {cdx:>06b} {cdx:03x} {is_fwd=!s:6s} {is_comp=!s:6s} {is_fake=!s:6s}")

		if kmer_size<7:
			ids[cdx] = ids.get(cdx, []) + [i]

		if cdx >= i:
			num_regs += 1
			if i % 100_000 == 0:
				print(f" {i=:15,d} / {4**kmer_size:15,d} | {num_regs=:15,d}")
			kmer.pack_w(cdx)

	print(f"{num_regs=:15,d}")

	kmer.close_w()

	print("created")

	return num_regs, ids

def check(kmer, num_regs=None, debug=False):
	kmer.open_r(list_regs=kmer.kmer_size<7, kmer_size=kmer.kmer_size, num_regs=num_regs)

	#kmer.list(kmer_size=kmer_size, list_regs=kmer_size<7, debug=kmer_size<7)

	print("checking")
	for i, seq in enumerate(generate_kmer(kmer.kmer_size)):
		if i % 100_000 == 0:
			print(f" {i=:15,d} / {4**kmer_size:15,d}")
		pos, val = kmer.find_sequence_id(seq)
		#cds, cdx, is_fwd, is_comp, is_fake = kmer.rev_comp_4(seq)
		#qes                                = kmer.generate_sequence(cdx)
		#if debug or kmer_size<7: print(f"{i=} {seq=} {pos=} {val=} {cds=} {cdx=} {qes=}")
		#assert qes == cds, f'{qes=} == {cds=}'
	print("checked")

	kmer.close_r()

def main(kmer_size, debug=False):
	kmer = Kmer(kmer_size=kmer_size, debug=kmer_size<7 or debug)

	num_regs = None
	if kmer.exists:
		print("exists")
	else:
		num_regs, ids = create(kmer, debug=debug)

	if kmer_size < 7:
		for i, (cdx, idxs) in enumerate(sorted(ids.items())):
			print(f"{i:3d} {cdx:3d} {idxs}")

	check(kmer, num_regs=num_regs, debug=debug)


if __name__ == "__main__":
	kmer_size = int(sys.argv[1])
	main(kmer_size=kmer_size)






















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

def rev_comp(seq):
	return "".join([RC[c] for c in seq[::-1]])

