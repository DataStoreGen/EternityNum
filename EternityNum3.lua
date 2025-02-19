-- made by gojosatoru7304 or DataStore_Genius
local ModuleCleaner = require(script.Parent.ModuleCleaner)
local expl = 1e10
local ldown = math.log10(expl)
local msd = 100
local AllowOverflow = true
local SuffixLimit = "9e1E14"
local DefaultDigits = 2
local tau = 6.2831853071795864769252842
local EXPN1 = 0.36787944117144232159553
local OMEGA = 0.56714329040978387299997

local st = {}
local C = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7}
local Sets = {"k","M","B"}
local FirstOnes = {"", "U","D","T","Qd","Qn","Sx","Sp","Oc","No"}
local SecondOnes = {"", "De","Vt","Tg","qg","Qg","sg","Sg","Og","Ng"}
local ThirdOnes = {"", "Ce", "Du","Tr","Qa","Qi","Se","Si","Ot","Ni"}
local MultOnes = {
	"", "Mi","Mc","Na","Pi","Fm","At","Zp","Yc", "Xo", "Ve", "Me", 
	"Due", "Tre", "Te", "Pt", "He", "Hp", "Oct", "En", "Ic", "Mei", 
	"Dui", "Tri", "Teti", "Pti", "Hei", "Hp", "Oci", "Eni", "Tra","TeC",
	"MTc","DTc","TrTc","TeTc","PeTc","HTc","HpT","OcT","EnT","TetC","MTetc",
	"DTetc","TrTetc","TeTetc","PeTetc","HTetc","HpTetc","OcTetc","EnTetc","PcT",
	"MPcT","DPcT","TPCt","TePCt","PePCt","HePCt","HpPct","OcPct","EnPct","HCt",
	"MHcT","DHcT","THCt","TeHCt","PeHCt","HeHCt","HpHct","OcHct","EnHct","HpCt",
	"MHpcT","DHpcT","THpCt","TeHpCt","PeHpCt","HeHpCt","HpHpct","OcHpct","EnHpct",
	"OCt","MOcT","DOcT","TOCt","TeOCt","PeOCt","HeOCt","HpOct","OcOct","EnOct","Ent","MEnT",
	"DEnT","TEnt","TeEnt","PeEnt","HeEnt","HpEnt","OcEnt","EnEnt","Hect", "MeHect"}

function F_gam(n)
	if n > 171.6236 then return 1.8e308 end
	if (n > 0.5) then
		n -= 1
		local x = C[1]
		for i=1,7 do
			x+=C[i+1]/(n+1)
		end
		local t = n+7.5
		return x*t^(n+0.5-36)*math.exp(-t)*t^36*2.50662827463100050241576528
	end
	return 3.141592653589793238 / (math.sin(3.141592653589793238 * n) * F_gam(1 - n))
end

function f_lamb(z)
	local tol, w, wn = 1e-10, nil, nil
	if z > 1.79e308 then return z elseif z == 0 then return z elseif z == 1 then return OMEGA end
	if z < 10 then
		w = 0 else w = math.log(z)-math.log(math.log(z))
	end
	for i=1, 100 do
		wn = (z + math.exp(-w)+w*w)/(w+1)
		if math.abs(wn-w) < tol * math.abs(wn) then return wn else w = wn end
	end
	error('Failed to itterate at: f_lamb')
end

function Cnew(man, metaExp, exp)
	return ModuleCleaner.new():RegisterTable('CnewFunc', {man = man, metaExp = metaExp, exp=exp})
end

local ZERO = Cnew(0, 0, 0)
local ONE = Cnew(1, 0, 1)
local NaN = Cnew(1, -1, 1)
local Inf = Cnew(1, math.huge, 1)
function st.IsNaN(Value) : boolean
	return Value.man == NaN.man and Value.metaExp == NaN.metaExp and Value.exp == NaN.exp 
end

function st.IsInf(Value) : boolean
	return Value.metaExp == math.huge or Value.exp == math.huge
end

function st.IsZero(Value) : boolean
	return Value.man == 0 or (Value.exp == 0 and Value.metaExp == 0)
end

function st.toStr(man, exp)
	return man .. 'e' .. exp
end

function newMeta(man, meta, exp)
	local result = ModuleCleaner.new():RegisterTable('Cnew', Cnew(man, meta, exp))
	return result
end

function st.correct(value)
	if st.IsNaN(value) then return NaN elseif st.IsInf(value) then return Inf elseif st.IsZero(value) then return ZERO end
	local man, meta, exp = value.man, value.metaExp, value.exp
	if meta == 0 and exp < 0 then
		exp, man = -exp, -man
	end
	if meta == 0 and exp < 1e-10 then
		return newMeta(man, meta+1, math.log10(exp))
	end
	local abs, sign = math.abs(exp), math.sign(exp)
	if abs >= expl then
		return newMeta(man, meta+1, sign*math.log10(abs))
	end
	while abs < ldown and meta > 0 do
		meta -=1
		exp = meta == 0 and math.pow(10, exp) or sign * math.pow(10, abs)
		abs, sign = math.abs(exp), math.sign(exp)
	end
	if meta == 0 and exp < 0 then
		exp, man = -exp, -man
	elseif exp == 0 then
		man = 0
	end
	return newMeta(man, meta, exp)
end

function st.new(man, metaExp, exp)
	return st.correct({man=man,metaExp =metaExp, exp=exp})
end

function st.fromNumber(value)
	local num = {}
	num.man = math.sign(value)
	num.metaExp = 0
	num.exp = math.abs(value)
	return st.correct(num)
end

function st.fromSci(value)
	local split = value:split('e')
	local man, exp = split[1], split[2]
	man, exp = tonumber(man), tonumber(exp)
	local sign = math.sign(man)
	local over = math.floor(math.log10(man))
	if over > 0 then
		man/=10^over
		exp+=over
	end
	if exp == 0 then return st.new(math.sign(man), 0, man) end
	if man == 0 then return ZERO end

	if exp < 0 then
		if exp < -100 then return ZERO end
		local exp2 = math.log10(man)+exp
		return st.correct(st.new(sign, 1, exp2))
	end
	local exp2 = math.log10(man)+exp
	local meta = 1
	if exp2 > expl then
		exp2 = math.log10(exp2)
		meta += 1
	end
	return st.correct(st.new(sign, meta, exp2))
end

function st.default(value)
	local split = value:split(';')
	local man = math.sign(tonumber(split[1]))
	if man == 0 then man = 1 end
	local meta = math.abs(tonumber(split[1]))
	local exp = tonumber(split[2])
	return st.correct(st.new(man, meta, exp))
end

function st.fromString(value)
	if value:find('e') and not value:find(';') then
		return st.fromSci(value)
	elseif value:find(';') then
		return st.default(value)
	end
	if value == 'NaN' then return NaN elseif value == 'Inf' then return Inf elseif value == '' then return ZERO end
	return st.fromNumber(tonumber(value))
end

function st.toString(value)
	if st.IsNaN(value) then return 'NaN' elseif st.IsInf(value) then return 'Inf' end
	return value.metaExp .. ';' .. value.exp
end

function st.convert(value)
	if typeof(value) == 'number' then return st.fromNumber(value)
	elseif typeof(value) == 'string' then return st.fromString(value)
	elseif typeof(value) == 'table' then
		if #value == 2 then
			local Bnum = value[1] .. 'e' .. value[2]
			return st.fromSci(Bnum)
		elseif #value == 3 then
			return st.correct(st.new(value[1], value[2], value[3]))
		elseif value.man then
			return st.correct(st.new(value.man, value.metaExp, value.exp))
		end
	end
end

function st.toNumber(value)
	value = st.convert(value)
	if value.metaExp > 1 then
		if math.sign(value.exp) == -1 then
			return value.man * 0
		end
		return value.man * 1.8e308
	end
	if value.metaExp == 0 then return value.man * value.exp
	elseif value.metaExp == 1 then
		return value.man * 10^value.exp
	end
	return math.log10(-1)
end

function st.abs(value)
	value = st.convert(value)
	if value.man == 0 then return ZERO end
	return st.new(1, value.metaExp, value.exp)
end

function st.maxAbs(val1, val2)
	val1, val2 = st.convert(val1), st.convert(val2)
	if st.cmpAbs(val1, val2) < 0 then return st.toString(val2) end
	return st.toString(val1)
end

function st.cmpAbs(val1, val2)
	local la = nil
	local lb = nil
	if val1.exp > 0 then la = val1.metaExp else la =- val1.metaExp end
	if val2.exp > 0 then lb = val2.metaExp else lb =- val2.metaExp end
	if la > lb then return 1 elseif la < lb then return -1 end
	if val1.exp > val2.exp then return 1 elseif val1.exp < val2.exp then return -1 end
	return 0
end

function st.cmp(val1, val2)
	if val1.man > val2.man then return 1 elseif val1.man < val2.man then return -1 end
	return val1.man *  st.cmpAbs(val1, val2)
end

function st.le(val1, val2)
	return st.cmp(st.convert(val1), st.convert(val2)) == -1
end

function st.me(val1, val2)
	return st.cmp(st.convert(val1), st.convert(val2)) == 1
end

function st.eq(val1, val2)
	return st.cmp(st.convert(val1), st.convert(val2)) == 0
end

function st.leeq(val1, val2)
	return not (st.cmp(st.convert(val1), st.convert(val2)) == 1)
end

function st.meeq(val1, val2)
	return not (st.cmp(st.convert(val1), st.convert(val2)) == -1)
end

function st.recip(value)
	value = st.convert(value)
	if value.exp == 0 then return NaN end
	if value.metaExp == 0 then return st.new(value.man, 0, 1/value.exp) end
	return st.new(value.man, value.metaExp, -value.exp)
end

function baseLog(value, base)
	value, base = st.convert(value), st.convert(base)
	if value.man <= 0 or base.man <= 0 then return NaN elseif st.IsNaN(base) or st.IsNaN(value) then return NaN end
	if value.metaExp == 0 and base.metaExp == 0 then return st.new(value.man, 0, math.log(value.exp)/math.log(base.exp)) end
	return st.div(st.log10(value), st.log10(base))
end

function st.neg(val)
	val = st.convert(val)
	return st.new(-val.man, val.metaExp, val.exp)
end

function st.log(value, base)
	if base then return baseLog(value, base) end
	value = st.convert(value)
	if value.man <= 0 then return NaN end
	if value.metaExp == 0 then
		return st.toString(st.new(value.man, 0, math.log10(value.exp)))
	elseif value.metaExp == 1 then
		local scale = 2.302585092994046
		return st.toString(st.new(math.sign(value.exp), 0, math.abs(value.exp) * scale))
	elseif value.metaExp == 2 then
		local scale = 0.36221568869946325
		return st.toString(st.new(math.sign(value.exp), 1, math.abs(value.exp)+scale))
	end
	return st.toString(st.new(math.sign(value.exp), value.metaExp-1, math.abs(value.exp)))
end

function st.log10(value)
	value = st.convert(value)
	if value.man <= 0 then return NaN end
	if value.metaExp > 0 then return st.toString(st.new(math.sign(value.exp), value.metaExp-1, math.abs(value.exp))) end
	return st.toString(st.new(value.man, 0, math.log10(value.exp)))
end

function st.add(val1, val2)
	val1, val2 = st.convert(val1), st.convert(val2)
	if st.IsInf(val1) or st.IsInf(val2) then return Inf elseif st.IsZero(val1) then return val2 elseif st.IsZero(val2) then return val1 end
	if val1.man == -val2.man and val1.metaExp == val2.metaExp and val1.exp == val2.exp then return ZERO end
	local a, b
	if val1.metaExp >= 2 or val2.metaExp >= 2 then return st.toString(st.maxAbs(val1, val2)) end
	if st.cmpAbs(val1, val2) > 0 then a = val1; b = val2 else a = val2; b = val1 end
	if a.metaExp == 0 and b.metaExp == 0 then return st.toString(st.fromNumber(a.man*a.exp + b.man*b.exp)) end
	local la, lb = a.metaExp* math.sign(a.exp), b.metaExp * math.sign(b.exp)
	if la - lb >= 2 then return a end
	if la == 0 and lb == -1 then
		if math.abs(b.exp-math.log10(a.exp)) > msd then return a else
			local mag = 10^a.exp-math.log10(b.exp)
			local man = b.man + a.man*mag
			local result = st.new(math.sign(man), 1, math.log10(b.exp)+math.floor(math.log10(man)))
			return st.toString(result)
		end
	elseif la == 1 and lb == 0 then
		if math.abs(a.exp- math.log10(b.exp)) > msd then return a end
		local mag = 10^a.exp-math.log10(b.exp)
		local man = b.man + a.man*mag
		local result = st.new(math.sign(man), 1, math.log10(b.exp)+math.floor(math.log10(man)))
		return st.toString(result)
	end
	if math.abs(a.exp-b.exp) > msd then return a end
	local mag = 10^(a.exp - b.exp)
	local man = b.man+a.man*mag
	local result = st.new(math.sign(man), 1, math.log10(b.exp)+math.floor(math.log10(man)))
	return st.toString(result)
end

function st.mul(val1, val2)
	val1, val2 = st.convert(val1), st.convert(val2)
	if st.IsInf(val1) or st.IsInf(val2) then return Inf elseif st.IsZero(val1) or st.IsZero(val2) then return ZERO end
	if val1.metaExp == val2.metaExp and val1.exp == -val2.exp then return st.toString(st.new(val1.man * val2.man, 0, 1)) end
	local a, b
	if (val1.metaExp > val2.metaExp) or (val1.metaExp == val2.metaExp and math.abs(val1.exp) > math.abs(val2.exp)) then
		a = val1; b = val2
	else a = val2; b = val1
	end
	if a.metaExp == 0 and b.metaExp == 0 then
		return st.toString(st.fromNumber(a.man * b.man* a.exp * b.exp))
	elseif a.metaExp >= 3 or (a.metaExp - b.metaExp >= 2) then
		return st.toString( st.new(a.man * b.man, a.metaExp, a.exp))
	elseif a.metaExp == 1 and b.metaExp == 0 then
		return st.toString(st.new(a.man*b.man, 1, a.exp + math.log10(b.exp)))
	elseif a.metaExp == 1 and b.metaExp == 1 then
		return st.toString(st.new(a.man*b.man, 1, a.exp+b.exp))
	end
	if (a.metaExp == 2 and b.metaExp==1) or (a.metaExp ==2 and b.metaExp ==2) then
		local t = st.new(math.sign(b.exp), b.metaExp-1, math.abs(b.exp))
		local n = st.add(st.new(math.sign(a.exp), a.metaExp-1, math.abs(a.exp)), t)
		return st.toString(st.new(a.man * b.man, n.metaExp+1, n.man*n.exp))
	end
	return NaN
end

function st.div(val1, val2)
	return st.mul(st.convert(val1), st.recip(st.convert(val2)))
end

function st.sub(val1, val2)
	return st.add(st.convert(val1), st.neg(val2))
end

function st.abslog10(val)
	val = st.convert(val)
	if st.IsZero(val) then return NaN end
	if val.metaExp > 0 then return st.toString(st.new(math.sign(val.exp), val.metaExp-1, math.abs(val.exp))) end
	return st.toString(st.new(1, 0, math.log10(math.abs(val.exp))))
end

function st.pow10(val)
	val = st.convert(val)
	if st.IsInf(val) then return Inf end
	if val.metaExp == 0 then
		local n = 10^(val.man*val.exp)
		if n<1.8e308 and math.abs(n) > 0.1 then return st.toString(st.new(1,0, n))
		else
			if val.man == 0 then return ONE end
			val = st.toString(st.new(val.man, val.metaExp+1, math.log10(val.exp)))
		end
	elseif val.man > 0 and val.exp > 0 then return st.toString(st.new(val.man, val.metaExp+1, val.exp))
	elseif val.man < 0 and val.exp > 0 then return st.toString(st.new(-val.man, val.metaExp+1, -val.exp))
	end
	return ONE
end

function st.pow(val1, val2)
	val1, val2 = st.convert(val1), st.convert(val2)
	if st.IsZero(val1) then return ZERO end
	if val1.man == 1 and val1.metaExp == 0 and val1.exp == 1 then return ONE
	elseif st.IsZero(val2) then return ONE
	elseif val2.man == 1 and val2.metaExp == 0 and val2.exp == 1 then return val1
	end
	local calc = st.pow10(st.mul(st.abslog10(val1), val2))
	if val1.man == -1 and st.toNumber(val1) % 2 == 1 then return st.neg(calc)
	elseif val1.man == -1 and st.toNumber(val2) < 1e20 then
		local comp = st.fromNumber(math.cos(st.toNumber(val2)*math.pi))
		return st.mul(calc, comp)
	end
	return calc
end

function st.toScience(value)
	if value.metaExp > 2 then if AllowOverflow then return '' end return 'Inf'
	elseif value.metaExp == 2 and value.exp > 308 then return 'Inf'
	elseif st.IsZero(value) then return '0e0'
	end
	if value.metaExp == 0 then
		local man = (value.exp/10^math.floor(math.log10(value.exp))) * value.man
		return st.toStr(man, math.floor(math.log10(value.exp)))
	elseif value.metaExp == 1 then
		local man = (10^(value.exp-math.floor(value.exp)))*value.man
		return st.toStr(man, math.floor(value.exp))
	end
end

function st.toSci(value)
	local exp = math.floor(math.log10(math.abs(value)))
	local man = value /10^exp
	return man .. 'e' .. exp
end

function CutDig(value, digit)
	if digit < 0 then return value end
	return math.floor(value * 10^digit)/10^digit
end

function st.toSuffix(value, digit: number?)
	digit = digit or DefaultDigits
	local bnum = st.toScience(value):split('e')
	local man: number, exp: number = bnum[1], bnum[2]
	local SNum = math.fmod(exp, 3)
	exp = math.floor(exp/3)-1
	if exp <= -1 then return CutDig(bnum[1]*10^bnum[2], digit)
	elseif exp < 3 then return CutDig(man*10^SNum, digit) .. Sets[exp+1]
	end
	local txt = ''
	local function suffix1(n)
		local hund = math.floor(n/100)
		n = math.fmod(n, 100)
		local tens = math.floor(n/10)
		n = math.fmod(n, 10)
		local one = math.floor(n/1)
		txt = txt .. FirstOnes[one+1]
		txt = txt .. SecondOnes[tens+1]
		txt = txt .. ThirdOnes[hund+1]
	end
	local function suffix2(n)
		if n > 0 then
			n += 1
		elseif n > 1000 then
			n = math.fmod(n, 1000)
		end
		suffix1(n)
	end
	if exp < 1000 then
		suffix1(exp)
		return CutDig(man*10^SNum, digit) .. txt
	end
	for i = math.floor(math.log10(exp)/3), 0, -1 do
		if exp >= 10 ^(i*3) then
			suffix2(math.floor(exp/10^(i*3))-1)
			txt = txt .. MultOnes[i+1]
			exp = math.fmod(exp, 10^(i*3))
		end
	end
	return CutDig(man*10^SNum, digit) .. txt
end

function st.between(val1, x, y)
	val1, x, y = st.convert(val1), st.convert(x), st.convert(y)
	return st.me(val1, x) and st.le(val1, y)
end

function st.toLayer(value, digits: number?)
	value = st.convert(value)
	digits = digits or DefaultDigits
	if st.between(value, ZERO, ONE) then
		local den = st.short(st.div(ONE, value))
		return '1/' .. den
	end
	if value.man == 1 then
		if value.exp < 0 then
			return 'E(' .. value.metaExp .. '-)' .. CutDig(math.abs(value.exp), digits)
		end
		return 'E(' .. value.metaExp .. ')' .. CutDig(value.exp, digits)
	end
	if value.man == 0 then return 'E(0)0' end
	return st.toLayer(st.abs(value), digits)
end

function st.short(value, digits: number?)
	value = st.convert(value)
	if st.le(value, SuffixLimit) then return st.toSuffix(value, digits) end
	return st.toLayer(value, digits)
end

function st.root(val1, val2)
	return st.pow(st.convert(val1), st.recip(st.convert(val2)))
end

function st.sqrt(value)
	return st.root(st.convert(value), 2)
end

function st.shift(value, digits)
	value = st.convert(value)
	if value.metaExp > 1 then return value
	elseif digits > 20 then return value
	end
	local d = 10^(value.exp-math.floor(value.exp))
	d = math.floor(digits*10^digits)/10^digits
	value.exp = math.floor(value.exp) + math.log10(d)
	return st.toString(value)
end

function st.exp(value)
	value = st.convert(value)
	if value.metaExp == 0 then
		if value.exp <= 709.7 then
			return st.fromNumber(math.exp(value.man*value.exp))
		end
		return st.toString(st.new(1, 1, value.man*math.log10(math.exp(1))*value.exp))
	end
	return st.toString(st.new(1, value.metaExp + 1, value.metaExp == 1 and 
		value * (math.log10(math.exp(1)) + value.exp) or value.man + value.exp))
end

function st.gamma(value)
	value = st.convert(value)
	if st.leeq(value, ZERO) then return NaN elseif value.exp < 0 then return st.recip(value) end
	if value.metaExp == 0 then if st.le(value, {1, 0, 24}) then	return st.fromNumber(F_gam(value.man * value.exp)) end
		local t = value.exp - 1
		local l = 0.9189385332046727
		l = (l + ((t + 0.5) * math.log(t)))
		l = l - t
		local n2 = t * t
		local np = t
		local lm = 12 * np
		local adj = 1 / lm
		local l2 = l + adj
		if (l2 == l) then return st.exp(l) end
		l = l2
		np = np * n2
		lm = 360 * np
		adj = 1 / lm
		l2 = l - adj
		if l2 == l then return st.exp(l) end
		l = l2
		np = np * n2
		lm = 1260 * np
		local lt = 1 / lm
		l = l + lt
		np = np * n2
		lm = 1680 * np
		lt = 1 / lm
		l = l - lt
		return st.exp(l)
	elseif value.metaExp == 1 then 
		return st.exp(st.mul(value, st.sub(st.log(value), 1)))
	end
	return st.exp(value)
end

function st.lbencode(value)
	value = st.convert(value)
	if st.eq(value, 1) then return 1 elseif value.man == 0 then return 4e18 end
	local man: number, meta: number, exp: number = value.man, value.metaExp, value.exp
	local mode = nil
	if man == -1 then
		if meta > 999 then
			mode = (math.sign(exp) == 1) and 0 or 2
		else
			mode = (math.sign(exp) == 1) and 1 or 3
		end
	elseif man == 1 then
		if meta > 9999 then
			mode = (math.sign(exp) == 1) and 8 or 6
		else
			mode = (math.sign(exp) == 1) and 7 or 5
		end
	end
	local v = mode * 1e18
	local logMeta = math.log10(meta)
	local logExp = math.log10(math.abs(exp):: number)
	local scaler = 3.2440674117208e15
	if mode == 8 then
		v += (logMeta+(logExp/10)) * scaler
	elseif mode == 7 then
		v += meta * 1e14 + logExp * 1e13
	elseif mode == 6 then
		v += 1e18-(logMeta+(logExp/10)) * scaler
	elseif mode == 5 then
		v += meta * 1e14 + 1e14 - logExp * 1e13
	elseif mode == 3 or mode == 2 then
		local off = 1e18-((logMeta+(logExp/10)) * scaler)
		v += off
	elseif mode == 1 then
		local off = 1e18 - (logMeta + logExp) * 1e13
		v += off
	elseif mode == 0 then
		local off = (logMeta + (logExp/10)) * scaler
		off = 1e18 - off
		v += off
	end
	return v
end

function st.lbdecode(value)
	if value == 2e18 then return st.new(-1, 0, 1) end
	if value == 3e18 then return st.new(-1, 10000, -1) end
	if value == 1e18 then return st.new(-1, 0, -1) end
	if value == 6e18 then return ONE end
	if value == 7e18 then return st.new(1, 10000, 1) end
	if value == 5e18 then return st.new(1, 10000, -1) end
	if value == 1 then return ONE end
	local mode: number = math.floor(value:: number / 1e18)
	if mode == 4 then return ZERO end
	local function decode(offset, neg)
		local v: number = value - offset
		if neg then v = 1e18 - v end
		v = v / 3.2440674117208e15
		local meta = math.floor(v)
		local n = 10^(math.fmod(v, 1) * 10)
		return meta, n
	end
	local meta, n
	if mode == 0 then
		meta, n = decode(0, true)
		return st.toString(st.new(-1, meta, n))
	elseif mode == 8 then
		meta, n = decode(8e18, false)
		return st.toString(st.new(1, meta, n))
	elseif mode == 1 then
		local v = value - 1e18
		meta = math.floor(v / 1e14)
		n = 10^(math.fmod(v, 1e14) / 1e13)
		return st.toString(st.new(1, meta, n))
	elseif mode == 7 then
		local v = value - 7e18
		meta = math.floor(v / 1e14)
		n = 10^(math.fmod(v, 1e14) / 1e13)
		return st.toString(st.new(-1, meta, n))
	elseif mode == 2 then
		meta, n = decode(2e18, true)
		return st.toString(st.new(-1, meta, -n))
	elseif mode == 6 then
		meta, n = decode(6e18, false)
		return st.toString(st.new(1, meta, n))
	elseif mode == 5 then
		local v = value - 5e18
		meta = math.floor(v / 1e14)
		n = 10^((1e14 - math.fmod(v, 1e14)) / 1e13)
		return st.toString(st.new(1, meta, -n))
	elseif mode == 3 then
		local v = value - 3e18
		v = 1e18 - v
		meta = math.floor(v / 1e14)
		n = 10^((1e14 - math.fmod(v, 1e14)) / 1e13)
		return st.toString(st.new(-1, meta, -n))
	end
	return NaN
end

function st.fact(value)
	return st.gamma(st.add(st.convert(value), 1))
end

function st.max(c, b, r, k)
	local en = st
	local max = en.div(en.log(en.add(en.div(en.mul(c , en.sub(r , 1)) , en.mul(b , en.pow(r,k))) , 1)) , en.log(r))
	local cost =  en.mul(b , en.div(en.mul(en.pow(r,k) , en.sub(en.pow(r,max) , 1)), en.sub(r , 1)))
	local nextCost = en.mul(b, en.pow(r,max))
	return max, cost, nextCost
end

function st.percent(val1, val2)
	val1, val2 = st.convert(val1), st.convert(val2)
	local val3 = st.correct(val1)
	local val4 = st.correct(val2)
	local division = st.div(val3, val4)
	if st.leeq(division, st.fromNumber(0.00001)) then return '0%' end
	if st.meeq(division, st.fromNumber(1)) then
		return '100%'
	end
	local percent = st.mul(division, st.fromNumber(100))
	percent = st.short(percent)
	return string.format('%s%%', percent)
end

function st.random(min, max)
	local seed = math.random()
	local even = st.sub(max, min)
	even = st.mul(even, seed)
	return st.add(even, min)
end

function st.exporand(min, max)
	min, max = st.convert(min), st.convert(max)
	local s1, s2 = min.man, max.man
	min = st.mul(st.exp(st.abs(min)), s1)
	max = st.mul(st.exp(st.abs(max)), s2)
	return st.exp(st.random(min, max))
end

function st.Changed(value, func: (newValue: string, label: TextLabel|TextButton) ->(), label: TextLabel|TextButton)
	value.Changed:Connect(function(newValue)
		func(newValue, label)
	end)
end

ModuleCleaner.new():RegisterTable('EternityNum3', st)

return st