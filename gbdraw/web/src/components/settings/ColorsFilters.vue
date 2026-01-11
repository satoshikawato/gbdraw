<script setup lang="ts">
import { HelpTip } from '@/components/common'
import type { SpecificRule, PriorityRule, FilterMode } from '@/types'

defineProps<{
  paletteNames: string[]
  selectedPalette: string
  currentColors: Record<string, string>
  specificRules: SpecificRule[]
  newSpecRule: SpecificRule
  filterMode: FilterMode
  manualBlacklist: string
  manualWhitelist: string[]
  priorityRules: PriorityRule[]
  newPriorityRule: PriorityRule
  customFeature: { name: string; color: string }
}>()

const emit = defineEmits<{
  'update:selectedPalette': [value: string]
  'update:currentColors': [feat: string, color: string]
  'update:filterMode': [value: FilterMode]
  'update:manualBlacklist': [value: string]
  'update:newSpecRule': [key: keyof SpecificRule, value: string]
  'update:newPriorityRule': [key: keyof PriorityRule, value: string]
  'update:customFeature': [key: 'name' | 'color', value: string]
  'reset-colors': []
  'add-custom-feature': []
  'add-specific-rule': []
  'remove-specific-rule': [index: number]
  'add-priority-rule': []
  'remove-priority-rule': [index: number]
  'add-whitelist-item': [item: string]
  'remove-whitelist-item': [index: number]
}>()

const newWhitelistItem = defineModel<string>('newWhitelistItem', { default: '' })
</script>

<template>
  <div class="card bg-slate-50/50">
    <details>
      <summary>
        <i class="ph ph-palette"></i> Colors & Filters
        <HelpTip text="Customize colors and filter labels by keyword." />
      </summary>

      <div class="mt-4 space-y-6 pt-4 border-t border-slate-200">
        <!-- Default Colors -->
        <div>
          <div class="flex justify-between items-center mb-2">
            <h4 class="text-xs font-bold text-slate-500">
              DEFAULT COLORS (-d)
              <HelpTip text="Base colors for features. Acts as the -d override table." />
            </h4>
            <button @click="emit('reset-colors')" class="text-[10px] text-blue-600 hover:underline">
              Reset
            </button>
          </div>

          <select
            :value="selectedPalette"
            @change="emit('update:selectedPalette', ($event.target as HTMLSelectElement).value)"
            class="form-input mb-3"
          >
            <option v-for="p in paletteNames" :key="p" :value="p">{{ p }}</option>
          </select>

          <div class="grid grid-cols-5 gap-1.5 bg-white p-2 rounded border border-slate-200">
            <div
              v-for="(color, feat) in currentColors"
              :key="feat"
              class="flex flex-col items-center group relative cursor-pointer"
            >
              <input
                type="color"
                :value="color"
                @input="emit('update:currentColors', String(feat), ($event.target as HTMLInputElement).value)"
                class="w-8 h-8 p-0 border-0 rounded cursor-pointer transition-transform hover:scale-110"
              />
              <span class="text-[9px] text-slate-500 mt-0.5 truncate w-full text-center" :title="String(feat)">
                {{ feat }}
              </span>
            </div>
          </div>

          <!-- Add Custom Feature -->
          <div class="mt-2 p-2 bg-slate-50 rounded border border-slate-200 flex gap-2 items-end">
            <div class="flex-grow">
              <label class="text-[9px] font-bold text-slate-500 block mb-1">ADD FEATURE</label>
              <input
                type="text"
                :value="customFeature.name"
                @input="emit('update:customFeature', 'name', ($event.target as HTMLInputElement).value)"
                class="form-input text-xs py-1"
                placeholder="Feature name"
              />
            </div>
            <input
              type="color"
              :value="customFeature.color"
              @input="emit('update:customFeature', 'color', ($event.target as HTMLInputElement).value)"
              class="w-8 h-8 p-0 border-0 rounded"
            />
            <button @click="emit('add-custom-feature')" class="btn btn-primary btn-sm">
              <i class="ph ph-plus"></i>
            </button>
          </div>
        </div>

        <!-- Specific Rules -->
        <div>
          <h4 class="text-xs font-bold text-slate-500 mb-2">
            SPECIFIC RULES (-t)
            <HelpTip text="Color rules based on qualifier values. Use regex patterns." />
          </h4>

          <!-- Existing rules -->
          <div v-if="specificRules.length > 0" class="space-y-1 mb-2">
            <div
              v-for="(rule, idx) in specificRules"
              :key="idx"
              class="flex items-center gap-2 p-1.5 bg-white rounded border border-slate-200 text-xs"
            >
              <span class="font-bold text-slate-600">{{ rule.feat }}</span>
              <span class="text-slate-400">{{ rule.qual }}:</span>
              <span class="text-blue-600 font-mono truncate flex-1">{{ rule.val }}</span>
              <input type="color" :value="rule.color" disabled class="w-5 h-5 p-0 border-0 rounded" />
              <button @click="emit('remove-specific-rule', idx)" class="text-red-500 hover:text-red-700">
                <i class="ph ph-x"></i>
              </button>
            </div>
          </div>

          <!-- Add new rule -->
          <div class="grid grid-cols-12 gap-1.5 p-2 bg-slate-50 rounded border border-slate-200">
            <select
              :value="newSpecRule.feat"
              @change="emit('update:newSpecRule', 'feat', ($event.target as HTMLSelectElement).value)"
              class="form-input text-xs py-1 col-span-3"
            >
              <option value="CDS">CDS</option>
              <option value="tRNA">tRNA</option>
              <option value="rRNA">rRNA</option>
              <option value="ncRNA">ncRNA</option>
              <option value="misc_feature">misc_feature</option>
            </select>
            <select
              :value="newSpecRule.qual"
              @change="emit('update:newSpecRule', 'qual', ($event.target as HTMLSelectElement).value)"
              class="form-input text-xs py-1 col-span-3"
            >
              <option value="product">product</option>
              <option value="gene">gene</option>
              <option value="locus_tag">locus_tag</option>
              <option value="note">note</option>
            </select>
            <input
              type="text"
              :value="newSpecRule.val"
              @input="emit('update:newSpecRule', 'val', ($event.target as HTMLInputElement).value)"
              class="form-input text-xs py-1 col-span-4"
              placeholder="Regex pattern"
            />
            <input
              type="color"
              :value="newSpecRule.color"
              @input="emit('update:newSpecRule', 'color', ($event.target as HTMLInputElement).value)"
              class="w-full h-7 p-0 border-0 rounded col-span-1"
            />
            <button @click="emit('add-specific-rule')" class="btn btn-primary btn-sm col-span-1">
              <i class="ph ph-plus"></i>
            </button>
          </div>
        </div>

        <!-- Label Filtering -->
        <div>
          <h4 class="text-xs font-bold text-slate-500 mb-2">
            LABEL FILTERING
            <HelpTip text="Filter which labels are displayed on the diagram." />
          </h4>

          <div class="flex gap-2 mb-2">
            <label
              v-for="fm in ['None', 'Blacklist', 'Whitelist'] as const"
              :key="fm"
              class="flex items-center gap-1 text-xs cursor-pointer"
            >
              <input
                type="radio"
                :value="fm"
                :checked="filterMode === fm"
                @change="emit('update:filterMode', fm)"
                class="text-blue-600"
              />
              {{ fm }}
            </label>
          </div>

          <!-- Blacklist input -->
          <div v-if="filterMode === 'Blacklist'">
            <textarea
              :value="manualBlacklist"
              @input="emit('update:manualBlacklist', ($event.target as HTMLTextAreaElement).value)"
              class="form-input text-xs"
              rows="2"
              placeholder="Comma-separated keywords to hide"
            ></textarea>
          </div>

          <!-- Whitelist input -->
          <div v-if="filterMode === 'Whitelist'" class="space-y-2">
            <div class="flex flex-wrap gap-1">
              <span
                v-for="(item, idx) in manualWhitelist"
                :key="idx"
                class="bg-green-100 text-green-700 px-2 py-0.5 rounded text-xs flex items-center gap-1"
              >
                {{ item }}
                <button @click="emit('remove-whitelist-item', idx)" class="hover:text-red-600">
                  <i class="ph ph-x text-[10px]"></i>
                </button>
              </span>
            </div>
            <div class="flex gap-1">
              <input
                v-model="newWhitelistItem"
                type="text"
                class="form-input text-xs py-1 flex-1"
                placeholder="Add keyword to whitelist"
                @keyup.enter="emit('add-whitelist-item', newWhitelistItem); newWhitelistItem = ''"
              />
              <button
                @click="emit('add-whitelist-item', newWhitelistItem); newWhitelistItem = ''"
                class="btn btn-primary btn-sm"
              >
                <i class="ph ph-plus"></i>
              </button>
            </div>
          </div>
        </div>

        <!-- Qualifier Priority -->
        <div>
          <h4 class="text-xs font-bold text-slate-500 mb-2">
            QUALIFIER PRIORITY
            <HelpTip text="Define which qualifier to use for labels (in priority order)." />
          </h4>

          <!-- Existing rules -->
          <div v-if="priorityRules.length > 0" class="space-y-1 mb-2">
            <div
              v-for="(rule, idx) in priorityRules"
              :key="idx"
              class="flex items-center gap-2 p-1.5 bg-white rounded border border-slate-200 text-xs"
            >
              <span class="font-bold text-slate-600">{{ rule.feat }}:</span>
              <span class="text-blue-600 font-mono flex-1">{{ rule.order }}</span>
              <button @click="emit('remove-priority-rule', idx)" class="text-red-500 hover:text-red-700">
                <i class="ph ph-x"></i>
              </button>
            </div>
          </div>

          <!-- Add new rule -->
          <div class="flex gap-1.5 p-2 bg-slate-50 rounded border border-slate-200">
            <select
              :value="newPriorityRule.feat"
              @change="emit('update:newPriorityRule', 'feat', ($event.target as HTMLSelectElement).value)"
              class="form-input text-xs py-1 w-24"
            >
              <option value="CDS">CDS</option>
              <option value="tRNA">tRNA</option>
              <option value="rRNA">rRNA</option>
              <option value="gene">gene</option>
            </select>
            <input
              type="text"
              :value="newPriorityRule.order"
              @input="emit('update:newPriorityRule', 'order', ($event.target as HTMLInputElement).value)"
              class="form-input text-xs py-1 flex-1"
              placeholder="product,gene,locus_tag"
            />
            <button @click="emit('add-priority-rule')" class="btn btn-primary btn-sm">
              <i class="ph ph-plus"></i>
            </button>
          </div>
        </div>
      </div>
    </details>
  </div>
</template>
