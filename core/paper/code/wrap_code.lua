function Code(el)
  if FORMAT:match('latex') then
    local s=el.text:gsub('\\','\\textbackslash{}'):gsub('([_%%&#$])','\\%1')
    local parts={}
    for word in s:gmatch('%S+') do table.insert(parts,'\\seqsplit{'..word..'}') end
    return pandoc.RawInline('latex', '\\texttt{'..table.concat(parts,'\\ ')..'}')
  end
end
