# Minimal UCSC .2bit reader (versions 0 and 1). Returns sequence with N blocks and soft-mask (lowercase) applied.
import struct

class TwoBit:
	def __init__(self, path):
		self.handle = open(path, "rb")
		signature, version, count, _ = struct.unpack("<IIII", self.handle.read(16))
		if signature != 0x1A412743:
			raise ValueError("not a little-endian 2bit file")
		offset_format = "<Q" if version == 1 else "<I"
		offset_size = 8 if version == 1 else 4
		self.offsets = {}
		for _ in range(count):
			name_size = self.handle.read(1)[0]
			name = self.handle.read(name_size).decode()
			self.offsets[name] = struct.unpack(offset_format, self.handle.read(offset_size))[0]

	def _header(self, name):
		self.handle.seek(self.offsets[name])
		dna_size, n_count = struct.unpack("<II", self.handle.read(8))
		n_starts = struct.unpack(f"<{n_count}I", self.handle.read(4 * n_count))
		n_sizes = struct.unpack(f"<{n_count}I", self.handle.read(4 * n_count))
		mask_count = struct.unpack("<I", self.handle.read(4))[0]
		mask_starts = struct.unpack(f"<{mask_count}I", self.handle.read(4 * mask_count))
		mask_sizes = struct.unpack(f"<{mask_count}I", self.handle.read(4 * mask_count))
		self.handle.read(4)
		return dna_size, list(zip(n_starts, n_sizes)), list(zip(mask_starts, mask_sizes)), self.handle.tell()

	def length(self, name):
		self.handle.seek(self.offsets[name])
		return struct.unpack("<I", self.handle.read(4))[0]

	def fetch(self, name, start, end):
		dna_size, n_blocks, mask_blocks, dna_offset = self._header(name)
		start, end = max(0, start), min(end, dna_size)
		first_byte = start // 4
		last_byte = (end + 3) // 4
		self.handle.seek(dna_offset + first_byte)
		packed = self.handle.read(last_byte - first_byte)
		bases = []
		for byte in packed:
			for shift in (6, 4, 2, 0):
				bases.append("TCAG"[(byte >> shift) & 3])
		sequence = bases[start - first_byte * 4:start - first_byte * 4 + (end - start)]
		for block_start, block_size in n_blocks:
			for position in range(max(block_start, start), min(block_start + block_size, end)):
				sequence[position - start] = "N"
		for block_start, block_size in mask_blocks:
			for position in range(max(block_start, start), min(block_start + block_size, end)):
				sequence[position - start] = sequence[position - start].lower()
		return "".join(sequence)

if __name__ == "__main__":
	# Self-check on a hand-built 2bit: one sequence "ACGTNNacgt" (N block at 4-6, mask at 6-10)
	import os, tempfile
	sequence = "ACGTNNACGT"
	code = {"T": 0, "C": 1, "A": 2, "G": 3, "N": 0}
	packed = bytearray()
	padded = sequence + "T" * (-len(sequence) % 4)
	for k in range(0, len(padded), 4):
		value = 0
		for base in padded[k:k + 4]:
			value = (value << 2) | code[base]
		packed.append(value)
	name = b"chrT"
	header_size = 16 + 1 + len(name) + 4
	record = struct.pack("<II", len(sequence), 1) + struct.pack("<I", 4) + struct.pack("<I", 2)
	record += struct.pack("<I", 1) + struct.pack("<I", 6) + struct.pack("<I", 4) + struct.pack("<I", 0) + bytes(packed)
	data = struct.pack("<IIII", 0x1A412743, 0, 1, 0) + bytes([len(name)]) + name + struct.pack("<I", header_size) + record
	path = os.path.join(tempfile.gettempdir(), "twobit_selfcheck.2bit")
	open(path, "wb").write(data)
	reader = TwoBit(path)
	assert reader.fetch("chrT", 0, 10) == "ACGTNNacgt", reader.fetch("chrT", 0, 10)
	assert reader.fetch("chrT", 3, 8) == "TNNac"
	print("twobit self-check ok")
