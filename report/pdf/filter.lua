-- filter.lua: pandoc Lua filter for the report's HTML (then printed to PDF by Chromium).
--  * the paragraph right after the title is the byline, the next one the licence line;
--  * "**Table N**" / "**Figure N**" label paragraphs and the italic title after them are kept
--    with the exhibit that follows (classes exhibit-label, exhibit-title);
--  * "*Note.*" paragraphs get class exhibit-note;
--  * paragraphs under "References" get a hanging indent (class references);
--  * bare URLs become links;
--  * table column widths follow the length of their contents (pipe tables otherwise get equal
--    widths).

local function text_of(inlines) return pandoc.utils.stringify(inlines) end

local function is_label(blk)
  if blk.t ~= "Para" or #blk.content ~= 1 or blk.content[1].t ~= "Strong" then return false end
  local s = text_of(blk.content[1].content)
  return s:match("^Table [%w%.]+$") ~= nil or s:match("^Figure %d+$") ~= nil
end
local function is_title(blk)
  return blk.t == "Para" and #blk.content == 1 and blk.content[1].t == "Emph"
end
local function is_note(blk)
  return blk.t == "Para" and blk.content[1] and blk.content[1].t == "Emph" and
         text_of(blk.content[1].content) == "Note."
end
local function div(blk, class) return pandoc.Div({ blk }, pandoc.Attr("", { class })) end

function Pandoc(doc)
  local out, blocks = {}, doc.blocks
  local seen_title, after_title, in_refs = false, 0, false
  local i = 1
  while i <= #blocks do
    local b = blocks[i]
    if b.t == "Header" then
      if b.level == 1 then seen_title = true; after_title = 0 end
      in_refs = (b.level == 2 and text_of(b.content) == "References")
      table.insert(out, b)
    elseif seen_title and after_title < 2 and b.t == "Para" then
      after_title = after_title + 1
      table.insert(out, div(b, after_title == 1 and "byline" or "licence"))
    elseif is_label(b) then
      -- label, title, the table or figure, and its note stay together (class exhibit)
      local head = { div(b, "exhibit-label") }
      if blocks[i + 1] and is_title(blocks[i + 1]) then
        table.insert(head, div(blocks[i + 1], "exhibit-title")); i = i + 1
      end
      local group = { pandoc.Div(head, pandoc.Attr("", { "exhibit-head" })) }
      local nxt = blocks[i + 1]
      if nxt and (nxt.t == "Table" or (nxt.t == "Para" and nxt.content[1] and nxt.content[1].t == "Image")) then
        table.insert(group, nxt); i = i + 1
        if blocks[i + 1] and is_note(blocks[i + 1]) then
          table.insert(group, div(blocks[i + 1], "exhibit-note")); i = i + 1
        end
      end
      table.insert(out, pandoc.Div(group, pandoc.Attr("", { "exhibit" })))
    elseif is_note(b) then
      table.insert(out, div(b, "exhibit-note"))
    elseif in_refs and b.t == "Para" then
      table.insert(out, div(b, "references"))
    else
      table.insert(out, b)
    end
    if b.t ~= "Header" and b.t ~= "Para" then after_title = 2 end
    i = i + 1
  end
  doc.blocks = out
  return doc
end

function Str(el)
  local url = el.text:match("^(https?://[^%s]+[^%s%.,;:%)])")
  if url then
    local rest = el.text:sub(#url + 1)
    local link = pandoc.Link(url, url)
    if rest == "" then return link end
    return { link, pandoc.Str(rest) }
  end
end

-- Column widths: proportional to the mean text length of the column's cells, but at least the
-- longest word, so that short label columns stay narrow and text columns get the room.
function Table(tbl)
  local ncol = #tbl.colspecs
  local sum, cnt, longest = {}, {}, {}
  for j = 1, ncol do sum[j], cnt[j], longest[j] = 0, 0, 1 end
  local function scan(rows)
    for _, row in ipairs(rows) do
      for j, cell in ipairs(row.cells) do
        if j <= ncol then
          local s = pandoc.utils.stringify(cell.contents)
          sum[j] = sum[j] + (utf8.len(s) or #s); cnt[j] = cnt[j] + 1
          for w in s:gmatch("%S+") do longest[j] = math.max(longest[j], utf8.len(w) or #w) end
        end
      end
    end
  end
  scan(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do scan(body.body) end
  local w, total = {}, 0
  for j = 1, ncol do
    w[j] = math.max(sum[j] / math.max(cnt[j], 1), 1.3 * longest[j], 3)
    total = total + w[j]
  end
  for j = 1, ncol do tbl.colspecs[j][2] = w[j] / total end
  return tbl
end
