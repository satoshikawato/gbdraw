<script setup lang="ts">
import type { DiagramMode, InputType, FileState, LinearSeq } from '@/types'
import { FileUploader } from '@/components/common'

defineProps<{
  mode: DiagramMode
  cInputType: InputType
  lInputType: InputType
  files: FileState
  linearSeqs: LinearSeq[]
}>()

const emit = defineEmits<{
  'update:cInputType': [value: InputType]
  'update:lInputType': [value: InputType]
  'update:files': [key: keyof FileState, value: File | null]
  'update:linearSeq': [index: number, key: keyof LinearSeq, value: File | null]
  'add-linear-seq': []
  'remove-linear-seq': []
}>()
</script>

<template>
  <div class="card border-l-4 border-l-blue-500">
    <div class="card-header">
      <i class="ph ph-folder-open text-xl"></i> Input Genomes
    </div>

    <!-- Circular Mode -->
    <div v-if="mode === 'circular'" class="space-y-4">
      <div class="flex gap-4 text-sm font-bold bg-slate-50 p-2 rounded-lg inline-flex">
        <label class="flex items-center gap-2 cursor-pointer hover:text-blue-600">
          <input
            type="radio"
            :checked="cInputType === 'gb'"
            @change="emit('update:cInputType', 'gb')"
            class="text-blue-600 focus:ring-blue-500"
          />
          GenBank
        </label>
        <label class="flex items-center gap-2 cursor-pointer hover:text-blue-600">
          <input
            type="radio"
            :checked="cInputType === 'gff'"
            @change="emit('update:cInputType', 'gff')"
            class="text-blue-600 focus:ring-blue-500"
          />
          GFF3 + FASTA
        </label>
      </div>

      <div v-if="cInputType === 'gb'">
        <FileUploader
          label="GenBank File (.gb)"
          accept=".gb,.gbk,.txt"
          :modelValue="files.c_gb"
          @update:modelValue="emit('update:files', 'c_gb', $event)"
        />
      </div>
      <div v-else class="grid grid-cols-1 gap-3">
        <FileUploader
          label="GFF3 File"
          accept=".gff,.gff3"
          :modelValue="files.c_gff"
          @update:modelValue="emit('update:files', 'c_gff', $event)"
        />
        <FileUploader
          label="FASTA File"
          accept=".fasta,.fa,.fna"
          :modelValue="files.c_fasta"
          @update:modelValue="emit('update:files', 'c_fasta', $event)"
        />
      </div>
    </div>

    <!-- Linear Mode -->
    <div v-if="mode === 'linear'" class="space-y-4">
      <div class="flex gap-4 text-sm font-bold bg-slate-50 p-2 rounded-lg inline-flex mb-2">
        <label class="flex items-center gap-2 cursor-pointer hover:text-blue-600">
          <input
            type="radio"
            :checked="lInputType === 'gb'"
            @change="emit('update:lInputType', 'gb')"
            class="text-blue-600"
          />
          GenBank
        </label>
        <label class="flex items-center gap-2 cursor-pointer hover:text-blue-600">
          <input
            type="radio"
            :checked="lInputType === 'gff'"
            @change="emit('update:lInputType', 'gff')"
            class="text-blue-600"
          />
          GFF3 + FASTA
        </label>
      </div>

      <div
        v-for="(seq, idx) in linearSeqs"
        :key="idx"
        class="p-3 bg-slate-50 rounded-lg border border-slate-200 relative group text-sm transition-all hover:border-blue-300"
      >
        <div class="absolute top-2 right-2 text-[10px] font-bold text-white bg-slate-400 px-1.5 py-0.5 rounded-full">
          #{{ idx + 1 }}
        </div>

        <div v-if="lInputType === 'gb'" class="space-y-2">
          <FileUploader
            label="GenBank File"
            accept=".gb,.gbk"
            :modelValue="seq.gb"
            @update:modelValue="emit('update:linearSeq', idx, 'gb', $event)"
          />
        </div>
        <div v-else class="space-y-2 mb-2">
          <FileUploader
            label="GFF3"
            accept=".gff"
            :modelValue="seq.gff"
            @update:modelValue="emit('update:linearSeq', idx, 'gff', $event)"
          />
          <FileUploader
            label="FASTA"
            accept=".fasta"
            :modelValue="seq.fasta"
            @update:modelValue="emit('update:linearSeq', idx, 'fasta', $event)"
          />
        </div>

        <div
          v-if="idx < linearSeqs.length - 1"
          class="mt-2 pt-2 border-t border-slate-200 border-dashed"
        >
          <div class="text-[10px] text-slate-500 font-bold mb-1 flex items-center gap-1">
            <i class="ph ph-arrows-down-up"></i> Compare to next (BLAST)
          </div>
          <FileUploader
            label="BLAST TSV"
            accept=".txt,.tsv,.csv"
            :modelValue="seq.blast"
            @update:modelValue="emit('update:linearSeq', idx, 'blast', $event)"
          />
        </div>
      </div>

      <div class="flex gap-2">
        <button @click="emit('add-linear-seq')" class="btn btn-secondary text-xs w-full">
          <i class="ph ph-plus"></i> Add Seq
        </button>
        <button
          v-if="linearSeqs.length > 1"
          @click="emit('remove-linear-seq')"
          class="btn btn-danger text-xs w-full"
        >
          <i class="ph ph-minus"></i> Remove
        </button>
      </div>
    </div>
  </div>
</template>
