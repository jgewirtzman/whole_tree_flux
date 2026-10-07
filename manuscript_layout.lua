-- Select document parts before applying layout. All prose stays in manuscript.md.
local part = os.getenv('MANUSCRIPT_PART') or 'combined'
local function select_part(doc)
  local before, after = pandoc.List(), pandoc.List()
  local in_si = false
  local title, authors
  for _, block in ipairs(doc.blocks) do
    if not title and block.t == 'Header' then title = block.content end
    if not authors and block.t == 'Para' then authors = block end
    if block.t == 'Header' and block.level == 1 and
        pandoc.utils.stringify(block.content) == 'Supporting Information' then
      in_si = true
    end
    if in_si then after:insert(block) else before:insert(block) end
  end
  if #after == 0 then error('Missing Supporting Information heading') end
  if part == 'main' then doc.blocks = before
  elseif part == 'si' then
    after:insert(2, pandoc.Para(title))
    after:insert(3, authors)
    doc.blocks = after
  elseif part ~= 'combined' then error('Unknown MANUSCRIPT_PART: ' .. part) end
  return doc
end
local first_header = true
local function format_image(el)
  el.attributes['width'] = '6in'
  return el
end
local function format_header(el)
  if first_header then
    first_header = false
    return pandoc.Div({pandoc.Para(el.content)}, pandoc.Attr('', {}, {['custom-style']='Title'}))
  end
  if el.level == 1 and pandoc.utils.stringify(el.content) == 'Supporting Information' then
    return pandoc.Div({pandoc.Para(el.content)}, pandoc.Attr('', {}, {['custom-style']='SI Heading'}))
  elseif pandoc.utils.stringify(el.content) == 'Abstract' then
    return pandoc.Div({pandoc.Para(el.content)}, pandoc.Attr('', {}, {['custom-style']='Abstract Heading'}))
  end
  return el
end
local function format_para(el)
  if #el.content == 1 and el.content[1].t == 'Image' then
    return pandoc.Div({el}, pandoc.Attr('', {}, {['custom-style']='Figure'}))
  end
  local txt = pandoc.utils.stringify(el.content)
  if txt:match('^Table S%d%.') then
    local caption = pandoc.Div({el}, pandoc.Attr('', {}, {['custom-style']='Table Caption'}))
    if txt:match('^Table S1%.') then
      return pandoc.Div({el}, pandoc.Attr('', {}, {['custom-style']='Table Page Caption'}))
    end
    return caption
  end
end
return {{Pandoc = select_part}, {Image = format_image, Header = format_header,
  Para = format_para, HorizontalRule = function() return {} end}}
