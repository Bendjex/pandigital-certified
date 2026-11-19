# Revisión adversaria C10 — lista de problemas

Manuscrito: «Certified pandigital approximations of π, γ, δ and e» (`manuscript.md`).
Paquete: `pandigital_cert_v3_borrador.zip` (los 238 hashes de `SHA256SUMS` verifican).

## Qué he hecho y qué no

**Leído entero:** el manuscrito; `c5/ivl.py`, `qx.py`, `motorx.py`, `hitcheck.py`, `c5run.py`, `lvl.py`, `xlib.c`; `src/motora.c`, `motorc.c`, `exact.py`, `canon.py`, `motorb.py`, `verify.py`, `certify.py`, `resolve.py`; `c11/nm.py`; `demostraciones_e/e_F.py`; `analysis/`; `docs/motorx_changes.md`.

**Ejecutado aquí** (scripts en `pruebas_revisor/`):

| Prueba | Tamaño | Resultado |
|---|---|---|
| `ivl.fwd` encadenado frente a `lvl` y racionales exactos, árboles al azar de 2–8 hojas | 23 000 árboles; ~5 800 anclas, ~5 600 canales factorizados, ~80 offsets, 5 de tercer nivel | 0 contradicciones (un caso de tercer nivel comprobado a mano con mpmath) |
| Cobertura estructural: expresiones al azar con algún nodo fuera de rango frente a N̂(S) y a Q(S, N, ventana estrecha); conjuntos de 4, 5 y 6 dígitos, en orden y en cualquier orden (n = 6) | 10 934 valores de X(S); 1 262 con valor dentro de rango, 208 exactamente iguales a 1; 921 consultas | 0 sin cubrir |
| Piezas de `c5run.plan` frente a la separación entera (n = 6, t_N = 4) | 3 separaciones, 4,1 millones de hits | multiconjuntos idénticos |
| Tests del paquete: `t_anchor`, `t_factor`, `t_factor2`, `t_offset`, `t_ancpow`, `t_ancrat`, `t_fwd2`, `t_pre`, `t_pre2`, `t_level2`, `t_c2prop`, `t_pieces` (parte 1), `t4_random_expr` (n = 6 en orden, 30 ventanas) | — | 0 fallos |
| `exp_lo/exp_hi/ln_lo/ln_hi` de motor A frente a mpmath | 57 144 puntos | 0 fallos; λ₁ < log₂e < λ₂ |
| Tablas 1, 2, 3 y 5 frente a `certified.jsonl`, `e_7_1.jsonl`, `c11.jsonl`, con analizador y evaluador propios (dígitos, orden, rango de cada nodo, error a 15 cifras) | 59 + 59 casillas | coinciden, salvo lo indicado en §6 |
| `e_F.py` para 1–7 y 1–8; `e_records_bounds` recalculado aparte | — | reproducen las salidas guardadas |

**Reproducciones de extremo a extremo** (segunda tanda, con ficheros de resolución reconstruidos aquí con `build_res.sh`):

| Qué | Resultado |
|---|---|
| `res_9_1`, `res_7_0`, `res_8_0` y `res_9_0` reconstruidos (un hilo) | SHA-256 idénticos a los publicados en `resolution_files.sha256`; `res_6_0` con 4 hilos: mismo contenido ordenado |
| Teorema 1, n = 9 en orden (π, γ, δ) con `certify.py` | 3797 / 3224 / 4090 candidatos, 1 / 1 / 3 valores, mismos errores y fórmulas que `certified.jsonl` |
| Teorema 1, n = 7 en cualquier orden (π, γ, δ) | 299 / 309 / 315 candidatos, 1 / 2 / 1 valores, registros idénticos a los publicados |
| Teorema 2, n = 9 en orden, con las ventanas con y sin menos unario (`c5run.py`, por piezas) | 8 separaciones, 40 / 40 / 40 hits (8 en S₁ = 6789 y 32 en S₁ = 789, como en los registros), 0 sin decidir en ambas ventanas, 0 operaciones de motor B sin decidir; 5,4 min de CPU |
| Motor C (FLINT 3.0.1) frente a motor B, n = 9 en orden y n = 7 en cualquier orden, tres constantes; 256 bits, reglas optimizadas a 1024 bits con tablas a 128, y sin ellas | 0 caminos no vacíos descartados; a 1024 bits quedan exactamente los no vacíos |
| Teorema 2, δ, n = 8 en cualquier orden, con la ventana sin menos unario: la consulta que falta en el paquete (1.1). Pasos, registros y trazas en `pruebas_revisor/delta_n8_sin_menos/` | 127 separaciones de 127, 492 hits (los mismos, separación a separación, que los registros publicados con el radio estrecho), 0 sin decidir en ambas ventanas, 0 valores dentro de la ventana sin menos, 0 operaciones de motor B sin decidir; 20 min de CPU |

**No hecho:** no he re-ejecutado el teorema 1 para n = 8 y 9 en cualquier orden (sólo sus ficheros de resolución y los contadores de 1.4 y 4.1), ni el teorema 2 para esas casillas, salvo δ con n = 8 (1.1). No he repetido `e_F.py` para 1–9 ni `e_F_09.py`. Las páginas externas no son accesibles desde aquí (`erich-friedman.github.io` y `mathworld.wolfram.com` están fuera de los dominios de red permitidos): la columna «Previous» de la tabla 1, la tabla 7 y §3.3 quedan sin contrastar.

---

## 1. Teorema X estructural, filtros, hits → valores, variante sin menos unario

**1.1 — Hueco de demostración (afirmación sin respaldo en los datos). δ, n = 8, cualquier orden: motor X no se ejecutó con la ventana sin menos unario.**
- Dónde: leyenda de la tabla 6, última frase del teorema 2 (§6), §8, C.6; `results/c5/estado_c5.csv`.
- Por qué: en `estado_c5.csv` la fila de δ dice `separaciones_con_C11 = 0`, `radio_consulta = 2.2255e-12` (el radio de C7) y `estado_C11 = PENDIENTE (falta el optimo de C11)`. Los registros `registros_c5_sesion/n8/c5_8_0_p{0,1}.jsonl` tienen `Eq = 2.2254665888666076e-12` y ningún `E11` para δ, y `n9/radii.txt` da «n8 E11» sólo para π y γ. El óptimo sin menos de esa casilla es 1,20·10⁻¹¹ (tabla 5), 5,4 veces mayor que el radio consultado. La leyenda de la tabla 6 («each value is checked against both… 0 values left undecided for either variant») es falsa para esa casilla, y CONTEXTO y el TODO de la tabla 5 sólo dan por pendiente n = 9.
- Comprobado aquí: he lanzado la consulta que falta (`c5run.py 8 0 delta`, por piezas, con `E11_delta = 1.2001163303948682e-11`, el óptimo sin menos de `c11.jsonl`; radio de consulta 1,2001163316·10⁻¹¹; fichero de resolución reconstruido con el hash publicado). Resultado: 127 separaciones de 127, 492 hits, 0 sin decidir frente a las dos ventanas, 0 valores dentro de la ventana sin menos unario, 0 operaciones de motor B sin decidir; 20 min de CPU. Los hits por separación coinciden, las 127, con los de los registros publicados. Con el mismo código que cierra las demás casillas, la afirmación resulta cierta; lo que falta es su respaldo en el paquete.
- Propuesta: repetir esa ejecución en la máquina del artículo (pasos en `pruebas_revisor/delta_n8_sin_menos/pasos.txt`), añadir sus registros a `registros_c5_sesion/n8/`, regenerar `estado_c5.csv` (hoy da `E_C11` vacío y «PENDIENTE (C11 en el laboratorio)» para esa fila, aunque `c11.jsonl` ya trae el óptimo) y añadir «n8 E11 delta» a `radii.txt`. Hasta entonces, la leyenda de la tabla 6 no se sostiene con los datos entregados.

**1.2 — Sin problema encontrado: inducción y filtros.** No hay ningún filtro que descarte por pertenecer o no a V(S). Revisados: `finish_form` (`qx.py:1171–1216`: ventana, `possibly_in/out`, `c2_drop`), `resolve_vv`/`vv_hull` (`qx.py:1375–1518`), `refine_base1` (`qx.py:1718–1797`, sólo sobre sub-hits con `isv`), `out_parts`, `near1`, `zero_hits`, `motorx.N` y `resolve_sub` (`motorx.py:307–368`, `433–467`). Todos deciden sobre el valor del propio nodo o sobre una condición necesaria de la ventana. Las restricciones de segundo nivel (`pre2`, `pre_level2`, `level2_lin`) tienen signo estricto, y por eso `c2_drop` es correcto incluso con 1.3.

**1.3 — Redacción (enunciado de C.4.1). El signo S no es estricto.**
- Dónde: C.4.1, C.4.3(3), C.4.4(8), C.4.6(1); `ivl.py:270–272`, `348–351`, `1060–1061`.
- Por qué: las reglas «1 ± t» y «q ± t» ponen S ≠ 0 (o DS ≠ 0) con M = −∞ (o DM = −∞) cuando t puede ser 0 sin llevar la marca Z, por ejemplo t = |a − b| con a = b no detectado. Para u = 1 la afirmación «sign(ln u) = S» es falsa. `fwd(s+l; 1, [−∞, −5])` da S = 1, M = −∞. «If t = 0 … the flag O takes over» sólo vale si Z está marcado. No he encontrado ningún consumidor que lo explote.
- Propuesta: definir el campo como «S·ln u ≥ 0 y M ≤ ln|ln u| ≤ N, con ln 0 = −∞», y decir que ln u ≠ 0 sólo se deduce de M > −∞. Lo mismo para DS.

**1.4 — Hueco (defecto latente, condicionado). La caché de `members` no distingue la tabla refinada.**
- Dónde: `qx.py:944–969` (clave `(S, i)`), `motorx.py:235–270` (línea 247).
- Por qué: `refine_near1` llama a `self.members(S, e, i)` con los índices de la tabla antes de refinar; al refinar se quitan y añaden entradas, y la clave `(S, i)` pasa a designar otra entrada con la lista de miembros de la clase antigua. Afectaría a `resolve_sub`, `resolve_vv` y `hitcheck.tab`. Sólo ocurre si alguna clase tiene 0 en su enclosure de ln sin ser el 1 exacto.
- Evidencia: no ocurre en ninguna tabla V de las que usa motor X. Contadas con `check_amb.py`: n = 9 en orden hasta 5 dígitos (1 233 008 clases) y n = 7, 8 y 9 en cualquier orden hasta 4 dígitos (2,3, 4,9 y 9,3 millones de clases): 0 clases ambiguas. **Sin efecto en ninguna casilla, tampoco en la de 1.1.**
- Propuesta: clave por los campos de la entrada, o vaciar la caché tras refinar.

**1.5 — Redacción / paquete. Faltan ficheros citados.** `t_coinc.py` (C.6 y C.7(6), la prueba de cobertura estructural), `lab/c5cierre/` y `fase3.py` (README, CONTEXTO, `motorx_changes` 26) y `tests/t3_constants.py` (`motora.c:70`) no están en el zip.

**1.6 — Sin problema encontrado: hits → valores y C.6.** `hitcheck.py:66–111` evalúa todas las derivaciones (`uniq_all`) sobre los miembros exactos. Dos dependencias que conviene decir en el texto: (a) `lvl.Invalid` descarta en silencio (`hitcheck.py:133`), lo que es correcto porque los ceros exactos de `lvl` sólo usan «misma forma ⇒ mismo valor»; (b) n = 8 se cerró por piezas y la prueba de extremo a extremo de `plan` figura como «NO TERMINADO» en `resultados_tests.txt` (la mía, a escala n = 6, pasa). El argumento de C.6 y el hecho 4 son correctos; lo único que falla es 1.1.

**1.7 — Redacción (hueco menor). Lema 2.** El enunciado y la hipótesis de inducción sólo dan |val E′| = |val E|, pero el paso «−E₁′ tiene valor |val E₁| > 0» y el uso de E₂′ como exponente necesitan val E′ = val E. Además está enunciado para R(S) y C.2–C.3 lo usan para expresiones con nodos fuera de rango. La construcción sirve tal cual: enunciarlo con igualdad de valores y para toda expresión válida.

## 2. Reglas de motor X (C.4, `motorx_changes.md`)

No he encontrado ninguna desigualdad falsa ni ningún redondeo hacia dentro con efecto. Notas:

**2.1 — Redacción (rigor formal). C.4.4(3), `ivl.py:296`.** `qq = succ(fl(Q1)·fl(Q2))` supone exacta la conversión de int64 a double. Con Q1 = Q2 = 2⁵³ + 1 da 2¹⁰⁶ + 2⁵⁴, que es Q1·Q2 − 1. La cota queda hacia dentro en ~10⁻³² y la holgura de `ln_hi` (≥ 10⁻¹⁴) la absorbe, pero el argumento escrito no lo cubre. Usar `_lnq_raw`, que sí trata Q ≥ 2⁵³.

**2.2 — Redacción (base de confianza, C.7).** «numpy… only through operations that are correctly rounded» no es exacto: `IndexT` selecciona candidatos con `np.log` (`qx.py:398–401`, `500–507`) y depende de que su error sea menor que 2⁻⁴⁰ relativo; `ipow` usa `np.log2` (`ivl.py:982`). Añadirlo a la lista o usar `ln_lo/ln_hi`.

**2.3 — Error de cifra. Estadísticas de tests que no coinciden con el paquete** (`results/c5/resultados_tests.txt` y mis ejecuciones): `t_factor2` 294 / 2706 / 1458 en el texto frente a 300 / 2700 / 1348; `t_anchor` 6000 frente a 3000; `t_offset` 33 622 / 18 292 / 176 frente a 19 950 / 6059 / 111; `t_level2` 30 000 frente a 20 000; `t_unique` 300 frente a 200; `t_near1` 267 932 frente a 441 784; `t_ancpow` 707 exactos frente a 727.

**2.4 — Redacción.** C.7 numera «1, 3, 4…». En §6, «|(.6 − .7)^(−(8^9)) − .1^(−(2^27))|» no es una escritura con 1–9 (2^27 abrevia 2^(3·(4+5)), y el bloque izquierdo va primero).

## 3. Lema 15 y motor C

Por lectura, cada regla que descarta es rigurosa: ventanas dobles hacia fuera (`dwindows`), `full_may_hit` con la holgura de `filter_full`, hojas estrechas con la comparación más laxa, agrupación por envolvente, bloques. Los eventos sin decidir conservan el camino. Y en ejecución (tabla de arriba): 0 descartes indebidos en n = 9 en orden y n = 7 en cualquier orden, con y sin las reglas optimizadas.

**3.1 — Redacción (robustez).** `decide_line` copia la línea a un búfer de 8192 bytes y parte en 16 tokens con 6 campos (`motorc.c:780–781`, `1032`); una línea más larga se truncaría sin aviso. Los caminos reales tienen como mucho 7 tokens. Añadir una comprobación explícita.

**3.2 — Redacción (paquete). Guiones que no corren tal como se entregan.** `src/c7_verify.sh:8` hace `cd` a una ruta absoluta de otra sesión y llama a `./motorc_new`, que `build_linux.sh` no genera. `src/make_c7.sh` usa `../c5/coremath/pow.c`, `../c7_templates/` y `/mnt/user-data/…`, que no existen en el paquete. `build_res4.sh` necesita `./motora4` y ningún guion lo compila. El README remite a los dos primeros.

## 4. Apéndice A

**4.1 — Hueco de demostración con defecto latente. Lema 9: `code_for` y `finish_real` con cota superior +∞.**
- Dónde: `motora.c:81–88`, `271–286`; fila «Interval [lo, hi]» y punto «Code» de la prueba del lema 9.
- Por qué: si la cota superior desborda (succ de un `cr_pow` infinito en `pow_hull`, línea 350, o `exp_hi` en `finish_from_ln`, línea 309) y lo ≤ succ(10³⁰⁰), `code_for` calcula `frexp(∞)` y devuelve 1, y `finish_real` devuelve `D_OK`. Comprobado: `finish_real(1e300, 9.99e299, +inf)` devuelve una clase de código 1 con H = [0x1.7e43c8800759ap+996, 0x1.7e43c8800759ep+996], que no contiene ni lo. El lema afirma H(r, k) ⊇ [lo, hi]. Hace falta que y·ln x abarque de 690,8 a 709,8, es decir, una base cuyo ln tenga una anchura relativa ≥ 2,7 % (clase ordinaria a menos de 4·10⁻⁵ de 1, o anclada con μ poco preciso por cancelación) y un exponente grande.
- Evidencia: con un contador en `finish_real` sobre las tablas reales (hasta 5 dígitos), 0 casos en n ≤ 9 en orden (1,1 millones de llamadas), n = 6 (37,3 millones), n = 7 (144,1 millones) y n = 8 (402,9 millones) en cualquier orden. En n = 9 en cualquier orden, 0 casos en las tablas de hasta 5 dígitos (947,1 millones de llamadas; salida en `pruebas_revisor/out_probe_inf_9_0_5.txt`). No afecta, por tanto, a ninguna casilla: el defecto es latente y el lema 9 se cumple en todas las tablas que usa el artículo, pero su prueba no cubre el caso.
- Propuesta: añadir `if (!isfinite(hi)) return D_SING_W;` (con los recuentos de arriba, tablas y candidatos tienen que salir idénticos) y tratar el caso hi = +∞ en la prueba del lema 9.

**4.2 — Redacción (hueco menor). A.7, paso 2.** La prueba necesita |v* − c| ≤ ε y no lo dice; `certify.py:132` lo anuncia en un comentario pero no filtra ni comprueba. En los 59 certificados el óptimo está estrictamente dentro de `E_window` (comprobado). Añadir la comprobación y la frase.

**4.2 bis — Redacción (determinismo).** `certify.py` ordena las escrituras del óptimo sólo por longitud (`sorted(..., key=len)` sobre un conjunto). Con escrituras de igual longitud el orden cambia de una ejecución a otra: en δ, n = 9 en orden, mi registro coincide con el publicado salvo en el orden de las tres fórmulas. Ordenar por (longitud, texto).

**4.3 — Redacción.** A.0 y «Reproducibility» dicen que el error se da con 25 cifras significativas. Para δ las bolas tienen radio ≥ 10⁻³⁰ (por ejemplo `[1.377579677254041e-14 +/- 4.44e-30]`): 15–16 cifras en n = 8–9.

**4.4 — Hueco menor (método de §8). `nm.py:36–41`, `96–100`.** Con base negativa, la validez se decide con `is_int`, que sólo reconoce enteros de tipo S. Un exponente entero de tipo N se trata como no entero: (7^(2^.5))^(8^(3/6)) vale 2401 (`canon` lo demuestra) pero `is_int` da `None`, y (.14 − .9)^((7^(2^.5))^(8^(3/6))) = −6,8·10⁻²⁸⁷ es una expresión sin menos, válida y dentro de rango. En ese ejemplo la bandera se recupera por otra derivación con los mismos dígitos (7^((2^.5)·(8^(3/6))), exponente racional exacto), y no veo efecto plausible en ningún óptimo. Pero «OR-ed over all derivations» no demuestra completitud: o se prueba que siempre existe una derivación de tipo S, o se marca como no resuelto cuando la bola del exponente contiene un entero.

## 5. Apéndice B

El lema B2 y la reducción a 4 + 4 cubren todos los casos (SB(2) = 9·ln 10 y SB(3) = 10⁹·ln 8 recalculados; U de un dígito forzado; bloques ≤ 3 dígitos).

**5.1 — Hueco menor. El lema B1 sólo cubre N ≥ 1** y F incluye t = −1. Para N = −M el error es (e/2M)(1 + 11/(12M) + …), mayor que e/2M. Hace falta para pasar de «|N| máximo» a «error mínimo» y para «unique as a value of N», que es de |N|.

**5.2 — Redacción.** «Distinct writings by block values: 16» cuenta los cuatro signos de A y C. Con exponentes impares (1–9: 3⁸⁴ y 9⁴²) la mitad de las líneas de `eF_9.out`, leídas literalmente como (U + A^X)^(C^Y), valen (1 − 1/N)^N ≈ 1/e; sólo están en F con s o t compensados. Aclarar qué se cuenta.

**5.3 — Redacción (rigor de la criba).** `e_F.py:277` preselecciona pares con puntos medios en coma flotante y tolerancia 10⁻⁹ antes del test de solape. Es seguro porque los radios relativos son ~10⁻⁷⁵ (medido: 3·10⁻⁷⁵ sobre 451 122 entradas en 1–7), pero el script no lo comprueba. Añadir un `assert`.

## 6. Cifras del texto frente a los datos

**6.1 — Error de cifra.** e en orden, n = 8: el error es 2,33499·10⁻²⁴ (−log₁₀ = 23,6317156612, tabla B1), o sea 2,33·10⁻²⁴ y no 2,34 (tablas 1 y 2).

**6.2 — Error de cifra.** §1.1: 18 457 734 525 360 901 453 873 570 es el número de cifras de N = 3^(2^85). Con la definición de B.1, la fórmula de Sabey es correcta hasta …569 decimales (log₁₀(N/e) = …569,394). La cifra de Bamberger (…164) sí es ⌊log₁₀(N/e)⌋, así que el párrafo mezcla dos convenciones.

**6.3 — Redacción.** §6 y C.8, «Up to 4 digits: no new value appears»: con C.2, 9^(8^7) ∈ NE({7, 8, 9}). `measure.py` cuenta valores dentro de rango que sólo se alcanzan con intermedios extremos; decirlo así.

**6.4 — Redacción (reproducibilidad).** `analysis/fig_d1.py:13, 63` dibuja n = 9 en cualquier orden como pendiente (marcadores huecos); en el PDF salen llenos. El script del paquete no es el que generó la figura 1.

**6.5 — Redacción.** §8 dice que el radio de la pasada única (n = 8 y 9) es una cota de un encuentro en el medio. Para δ, n = 8 el radio usado (3,058787·10⁻¹¹) es exactamente la ventana C7 de δ con n = 7. No afecta: el óptimo hallado (1,20·10⁻¹¹) queda dentro.

Todo lo demás cuadra: 21 filas en la tabla 1 y 20 casillas mejoradas; 22 de 59 sin cambio y, con n ≥ 7, sólo δ en orden y γ en cualquier orden con n = 8; tabla 4 (137 445 líneas, 18 371 / 30 904 / 22 911 distintas); tabla 6 (hits y CPU); 18 451 / 31 016 / 23 230; récords de e y tabla B1; ln ln N de B3 y B2; errores de §1.1; Aizenbud; 2^(3^7)/9^(5^4); el desplazamiento de δ.

## 7. §9.2

Está presentado como heurístico, y las cifras se reproducen (39 de 42, 0,82, los tres casos fuera, predicciones para n = 9). Afirmaciones que se pasan:

**7.1 — Redacción.** La media geométrica 0,82 se lee como «cerca de 1», pero una exponencial de media 1 tiene media geométrica e^(−γ) = 0,561. Los óptimos son en media geométrica 1,46 veces peores que lo que predice el propio modelo (z ≈ 1,9 si las casillas fueran independientes, y no lo son).

**7.2 — Redacción.** «no sign of further structure»: 15 de las 42 casillas usan Dₙ extrapolado; los tres casos fuera de la banda están entre las 27 medidas y los tres por arriba. Con eso y 7.1, la frase es más fuerte que la evidencia.

**7.3 — Redacción.** «the expected error is proportional to Δ = 0.05·c» sólo vale porque Dₙ es casi igual para las tres constantes (4,23, 4,18 y 3,33 millones con n = 6): es una observación empírica sobre la densidad, y hay que decirla.

**7.4 — Redacción.** «computed Dₙ exactly, in double precision» se contradice. «Nₙ, the exact number of expressions of Section 2.1»: la gramática admite −−E, luego hay infinitas; Nₙ cuenta un signo opcional por nodo (N₁ = 4).

## 8. Bibliografía y atribuciones

**8.1 — Hueco de fundamento. δ.** La pertenencia de δ a la bola de radio 10⁻³⁰ se apoya en OEIS y Briggs, que son cálculos numéricos sin cota de error demostrada. Es la única entrada de la base de confianza que no es una enclosure rigurosa. Conviene citar una cota asistida por ordenador con error demostrado (creo que la hay con cientos de cifras, de Burbanks, Osbaldestin y Thurlby; no he podido verificar la referencia) o declarar los resultados de δ condicionados a esos dígitos.

**8.2 — Redacción.** [5] y [14] se citan de segunda mano («as reported in [3]»); la columna de *Technology Review* es *Puzzle Corner* (ed. A. Gottlieb). §9.1, «it is expected that it does not», va sin fuente.

Las referencias [6]–[10], [17] y [18] coinciden con lo que recuerdo (sin comprobación en línea). No he podido contrastar [1]–[4], [11], [15] ni [16].

---

## Veredicto

Nada de lo encontrado invalida el teorema 1 ni el teorema 2 tal como están enunciados (con sus marcas TODO); los defectos latentes de 1.4 y 4.1 no se dan en ninguna de las tablas que usan las casillas (0 casos en los recuentos completos).
La última frase del teorema 2 (variante sin menos unario) para δ, n = 8, cualquier orden (1.1) es cierta según la consulta que he ejecutado aquí, pero el paquete no la respalda: falta esa ejecución en los datos y `estado_c5.csv` la da por pendiente.
