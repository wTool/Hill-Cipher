/** ----------------------------------

Name    : 希尔密码 - Hill Cipher
Version : 1.0

----------------------------------
*/
function MatrixSolver() {
	this.c_valid = 0;
	this.c_nonsquare = -1;
	this.c_singular = -2;
	this.c_wrongdimensions = -3;

	function get(matrixarray, i, j) {
		var row = matrixarray[i];
		return row[j];
	}

	function columns(matrixarray) {
		var row = matrixarray[0];
		return row.length;
	}

	function rows(matrixarray) {
		return matrixarray.length;
	}

	function show(matrixarray) {
		var s = "";
		for (var i = 0; i < rows(matrixarray); ++i) {
			for (var j = 0; j < columns(matrixarray); ++j) {
				var row = matrixarray[i];
				s += " " + row[j];
			}
			s += "\n\r";
		}
		alert(s);
	}

	function minor(matrixarray, i, j) {
		var m = [];
		for (var k = 0; k < matrixarray.length; ++k) {
			if (k == i) continue;
			var row = matrixarray[k];
			var newrow = [];
			for (var l = 0; l < row.length; ++l) {
				if (l == j) continue;
				newrow.push(row[l]);
			}
			m.push(newrow);
		}
		return m;
	}

	this.calcMultiplication = function(matrixA, matrixB) {
		var result = [];
		var columnsA = columns(matrixA);
		var rowsA = rows(matrixA);
		var columnsB = columns(matrixB);
		var rowsB = rows(matrixB);
		if (columnsA != rowsB) throw this.c_wrongdimensions;


		for (var i = 0; i < rowsA; ++i) {
			var newrow = [];
			for (var j = 0; j < columnsB; ++j) {
				var value = 0;
				for (var k = 0; k < rowsB; ++k) {
					value += get(matrixA, i, k) * get(matrixB, k, j);
				}
				newrow.push(value);
			}
			result.push(newrow);
		}
		return result;
	}

	this.calcTransponent = function(matrixarray) {
		var transponent = [];
		var cols = columns(matrixarray);
		for (var i = 0; i < cols; ++i) {
			var newrow = [];
			var rowscount = rows(matrixarray);
			for (var j = 0; j < rowscount; ++j) {
				newrow[newrow.length] = get(matrixarray, j, i);
			}
			transponent[transponent.length] = newrow;
		}
		return transponent;
	}

	this.calcScalar = function(matrixarray, scalarValue) {
		var scalar = [];
		var cols = columns(matrixarray);
		for (var i = 0; i < cols; ++i) {
			var newrow = [];
			var rowscount = rows(matrixarray);
			for (var j = 0; j < rowscount; ++j) {
				newrow.push(get(matrixarray, i, j) * scalarValue);
			}
			scalar.push(newrow);
		}
		return scalar;
	}

	this.calcDeterminant = function(matrixarray) {
		var columnsA = columns(matrixarray);
		var rowsA = rows(matrixarray);
		if (columnsA != rowsA) throw this.c_nonsquare;

		if (matrixarray.length == 1) return get(matrixarray, 0, 0);
		else if (matrixarray.length == 2) return get(matrixarray, 0, 0) * get(matrixarray, 1, 1) - get(matrixarray, 0, 1) * get(matrixarray, 1, 0);
		else {
			var det = 0;
			var cols = columns(matrixarray);
			for (var i = 0; i < cols; ++i) {
				det += Math.pow(-1, i) * get(matrixarray, 0, i) * this.calcDeterminant(minor(matrixarray, 0, i));
			}
			return det;
		}
	}

	this.calcInverse = function(matrixarray) {
		var columnsA = columns(matrixarray);
		var rowsA = rows(matrixarray);
		if (columnsA != rowsA) throw this.c_nonsquare;
		var detA = this.calcDeterminant(matrixarray);
		if (detA == 0) throw this.c_singular;
		var minorsmatrix = [];
		for (var i = 0; i < rowsA; ++i) {
			var newrow = [];
			for (var j = 0; j < columnsA; ++j) {
				var val = this.calcDeterminant(minor(matrixarray, i, j));
				val = val * Math.pow(-1, i + j);
				newrow.push(val);
			}
			minorsmatrix.push(newrow);
		}
		var transponentminors = this.calcTransponent(minorsmatrix);
		var scalarresult = this.calcScalar(transponentminors, 1 / detA);
		return scalarresult;
	}

	function norm(ta, m) {
		if (ta > m) ta = ta % m;
		if (ta < 0) ta = ta + m * (Math.floor(Math.abs(ta) / m) + 1);
		return ta;
	}

	this.calcScalarMod = function(matrixarray, scalarValue, modulus) {
		var scalar = [];
		var cols = columns(matrixarray);
		for (var i = 0; i < cols; ++i) {
			var newrow = [];
			var rowscount = rows(matrixarray);
			for (var j = 0; j < rowscount; ++j) {
				newrow.push(norm(get(matrixarray, i, j) * scalarValue, modulus));
			}
			scalar.push(newrow);
		}
		return scalar;
	}

	this.calcDeterminantMod = function(matrixarray, modulus) {
		var columnsA = columns(matrixarray);
		var rowsA = rows(matrixarray);
		if (columnsA != rowsA) throw this.c_nonsquare;

		//show(matrixarray);

		if (matrixarray.length == 1) return norm(get(matrixarray, 0, 0), modulus);
		else if (matrixarray.length == 2) {
			return norm(get(matrixarray, 0, 0) * get(matrixarray, 1, 1) - get(matrixarray, 0, 1) * get(matrixarray, 1, 0), modulus);
		} else {
			var det = 0;
			var cols = columns(matrixarray);
			for (var i = 0; i < cols; ++i) {
				det += norm(Math.pow(-1, i) * get(matrixarray, 0, i) * this.calcDeterminantMod(minor(matrixarray, 0, i), modulus), modulus);
			}
			return norm(det, modulus);
		}
	}

	this.calcInverseMod = function(matrixarray, modulus) {
		var columnsA = columns(matrixarray);
		var rowsA = rows(matrixarray);
		if (columnsA != rowsA) throw this.c_nonsquare;
		var detA = this.calcDeterminantMod(matrixarray, modulus);
		if (detA == 0) throw this.c_singular;
		var minorsmatrix = [];
		for (var i = 0; i < rowsA; ++i) {
			var newrow = [];
			for (var j = 0; j < columnsA; ++j) {
				var val = this.calcDeterminantMod(minor(matrixarray, i, j), modulus);
				val = norm(val * Math.pow(-1, i + j), modulus);
				newrow.push(val);
			}
			minorsmatrix.push(newrow);
		}
		var transponentminors = this.calcTransponent(minorsmatrix);

		var result = Calculate3312({
			"a": detA,
			"m": modulus
		});
		var scalarresult = this.calcScalarMod(transponentminors, result, modulus);
		return scalarresult;
	}
}

function Calculate3312(params) {
	
	var a = params["a"] === undefined ? 3 : params["a"];
	var m = params["m"] === undefined ? 26 : params["m"];
	var errors = {
		"e": "There is no modular multiplicative inverse for this integer"
	};
	
	var ta = a;
	if (ta > m) ta = ta % m;
	if (ta < 0) ta = ta + m * (Math.floor(Math.abs(ta) / m) + 1);

	var result = Calculate3299({
		"first": PCR.integer(ta.valueOf()),
		"second": PCR.integer(m.valueOf())
	});
	var res = Number(result.res.toString());
	if (res != 1) throw {
		"source": "a",
		"message": errors["e"]
	};
	else {
		var c2 = Number(result.coef2.toString());
		var tx = (c2 % m + m) % m;
	}

	return tx;
};

function Calculate3299(params) {
	
	var first = params["first"] === undefined ? "180" : params["first"];
	var second = params["second"] === undefined ? "150" : params["second"];
	
	var euklid = {
		gcd: PCR.integer(1),
		x: PCR.integer(0),
		y: PCR.integer(0)
	};

	function gcd(a, b, holder) {
		if (a.eq(0)) {
			holder.x = PCR.integer(0);
			holder.y = PCR.integer(1);
			return b;
		}
		var q = b.divmod(a);
		var d = gcd(q.r, a, holder);
		var tx = holder.x;
		holder.x = holder.y.sub(q.q.mul(holder.x));
		holder.y = tx;
		return d;
	}

	euklid.gcd = gcd(first.gt(second) ? second : first, first.gt(second) ? first : second, euklid);
	
	return {
		"res" : euklid.gcd,
		"coef1" : euklid.y,
		"coef2" : euklid.x,
	};
	
};

function calculateDisplay( params ){
	
	var alphabet = params["alphabet"] === undefined ? "abcdefghijklmnopqrstuvwxyz ,." : params["alphabet"];
	var key = params["key"] === undefined ? "best hill cipher" : params["key"];
	var text = params["text"] === undefined ? "the truth is out there" : params["text"];
	var method = params["method"] === undefined ? "0" : params["method"];
	var errors = {
		"key_invalid": "Wrong key. Key length should be square of integer.",
		"key_wrong": "Wrong key. Key should contain only alphabet symbols",
		"text_wrong": "Wrong text. Text should contain only alphabet symbols",
		"NonSquare": "Wrong key. Key matrix should be square",
		"Singular": "Wrong key. Key has degenerate matrix",
		"key_bad": "Wrong key. Key matrix determinant does not have modular multiplicative inverse",
		"method": {
			"c": "Encrypt",
			"d": "Decrypt"
		}
	};
	
	if (!(Math.sqrt(key.length) % 1 === 0)) throw {
		"source": "key",
		"message": errors["key_invalid"]
	};
	var patt = new RegExp("[^" + alphabet + "]");
	if (patt.test(key)) throw {
		"source": "key",
		"message": errors["key_wrong"]
	};
	if (patt.test(text)) throw {
		"source": "text",
		"message": errors["text_wrong"]
	};
	
	var alphabet_digital = [];
	for (var i = 0; i < alphabet.length; ++i) {
		alphabet_digital[alphabet[i]] = i;
	}

	var key_size = Math.sqrt(key.length);
	var key_digital = new Array();
	for (var i = 0; i < key_size; ++i) {
		var inner = new Array();
		for (var j = 0; j < key_size; ++j) {
			inner.push(alphabet_digital[key[i * 3 + j]]);
		}
		key_digital.push(inner);
	}

	var solver = new MatrixSolver();

	var detm = solver.calcDeterminantMod(key_digital, alphabet.length);
	
	try {
		Calculate3312({
			"a": detm,
			"m": alphabet.length
		});
	} catch (e) {
		throw {
			"source": "key",
			"message": errors["key_bad"]
		};
	}
	
	var transformed_text = "";

	//showme(key_digital);
	if (method == "d") key_digital = solver.calcInverseMod(key_digital, alphabet.length);
	//showme(key_digital);
	var vector = new Array();
	while (vector.length < key_size) vector.push(new Array());

	for (var t = 0; t < text.length; ++t) {
		if (t > 0 && t % key_size == 0) {
			transformed_text += transform(vector, key_digital);
			vector = new Array();
			while (vector.length < key_size) vector.push(new Array());
		}
		vector[t % key_size].push(alphabet_digital[text[t]]);
	}

	//padding
	if (vector[0].length > 0) {
		for (var i = 0; i < key_size; ++i)
		if (vector[i].length == 0) vector[i].push(alphabet_digital[text[text.length - 1]])
		transformed_text += transform(vector, key_digital);
	}
	
	function showme(record) {
		var property, propCollection = "";

		for (property in record) {
			propCollection += (property + ": " + record[property] + "\n");
		}
		alert(propCollection);
	}

	function transform(part, key_matrix) {
		try {
			var res = solver.calcMultiplication(key_matrix, part);
			var s = "";
			for (var i = 0; i < res.length; ++i)
			s += alphabet[res[i] % alphabet.length];
			return s;
		} catch (err) {
			if (err == solver.c_nonsquare) throw {
				"source": "matrixA",
				"message": errors["NonSquare"]
			};
			else if (err == solver.c_singular) throw {
				"source": "matrixA",
				"message": errors["Singular"]
			};
			else alert(err);
		}

	}
	return transformed_text;
}

function GenRandKey(){
    var keychars = "abcdefghijklmnopqrstuvwxyz";
    var chars = keychars.split("");
    var ret="";
    for(i=0; i<16; i++){
        index = Math.floor(chars.length*Math.random());
        ret += chars[index];
        chars.splice(index,1);
    } 
	return ret;
}
