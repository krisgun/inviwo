/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2012-2015 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************************/

#ifndef IVW_DESERIALIZER_H
#define IVW_DESERIALIZER_H



#include <inviwo/core/io/serialization/serializebase.h>
#include <inviwo/core/util/exception.h>
#include <inviwo/core/util/stdextensions.h>
#include <inviwo/core/util/stringconversion.h>
#include <inviwo/core/io/serialization/deserializationerrorhandler.h>
#include <inviwo/core/io/serialization/nodedebugger.h>
#include <type_traits>
#include <list>

namespace inviwo {

class Serializable;
class VersionConverter;
class InviwoApplication;
class FactoryBase;
template <typename T, typename K>
class ContainerWrapper;


class IVW_CORE_API Deserializer : public SerializeBase {
public:
    /**
     * \brief Deserializer constructor
     *
     * @param fileName path to file that is to be deserialized.
     * @param allowReference flag to manage references to avoid multiple object creation.
     */
    Deserializer(InviwoApplication* app, std::string fileName, bool allowReference = true);
    /**
     * \brief Deserializes content from the stream using path to calculate relative paths to data.
     *
     * @param stream Stream with content that is to be deserialized.
     * @param path A path that will be used to decode the location of data during
     * deserialization.
     * @param bool allowReference flag to manage references to avoid multiple object creation.
     */
    Deserializer(InviwoApplication* app, std::istream& stream, const std::string& path,
                 bool allowReference = true);

    void pushErrorHandler(BaseDeserializationErrorHandler*);
    BaseDeserializationErrorHandler* popErrorHandler();

    virtual ~Deserializer();

    // std containers
    /**
     * \brief Deserialize a vector
     *
     * Deserialize the vector that has pre-allocated objects of type T or allocated by deserializer.
     * A vector is identified by key and vector items are identified by itemKey
     *
     * eg. xml tree with key=Properties and itemKey=Property
     *
     *     <Properties>
     *          <Property identifier="enableMIP" displayName="MIP">
     *              <value content="0" />
     *          </Property>
     *          <Property identifier="enableShading" displayName="Shading">
     *              <value content="0" />
     *          </Property>
     *     <Properties>
     *
     * @param key vector key.
     * @param sVector vector to be deserialized.
     * @param itemKey vector item key
     */
    template <typename T>
    void deserialize(const std::string& key, std::vector<T*>& sVector, const std::string& itemKey);

    template <typename T, typename C>
    void deserialize(const std::string& key, std::vector<T*>& sVector, const std::string& itemKey,
                     C identifier);

    template <typename T>
    void deserialize(const std::string& key, std::vector<T>& sVector, const std::string& itemKey);

    template <typename T>
    void deserialize(const std::string& key, std::vector<std::unique_ptr<T>>& vector,
                     const std::string& itemKey);

    template <typename T>
    void deserialize(const std::string& key, std::list<T>& sContainer, const std::string& itemKeyr);
    
    /**
     * \brief  Deserialize a map
     *
     * Deserialize a map, which can have
     * keys of type K,
     * values of type V* (pointers)
     * and compare function C ( optional if
     * K primitive type, i.e., std::string, int, etc.,)
     * eg., std::map<std::string, Property*>
     *
     * eg. xml tree
     *
     *     <Properties>
     *          <Property identifier="enableMIP" displayName="MIP">
     *              <value content="0" />
     *          </Property>
     *          <Property identifier="enableShading" displayName="Shading">
     *              <value content="0" />
     *          </Property>
     *     <Properties>
     *
     * In the above xml tree,
     *
     * key                   = "Properties"
     * itemKey               = "Property"
     * param comparisionAttribute  = "identifier"
     * param sMap["enableMIP"]     = address of a property
     *         sMap["enableShading"] = address of a property
     *         where, "enableMIP" & "enableShading" are keys.
     *         address of a property is a value
     *
     * Note: If children has attribute "type", then comparisionAttribute becomes meaningless.
     *       Because deserializer always allocates a new instance of type using registered
     *       factories. eg.,
     *           <Processor type="EntryExitPoints" identifier="EntryExitPoints" reference="ref2" />
     *
     * @param key Map key or parent node of itemKey.
     * @param sMap  map to be deserialized - source / input map.
     * @param itemKey map item key of children nodes.
     * @param comparisionAttribute forced comparison attribute.
     */
    template <typename K, typename V, typename C, typename A>
    void deserialize(const std::string& key, std::map<K, V, C, A>& sMap, const std::string& itemKey,
                     const std::string& comparisionAttribute = SerializeConstants::KeyAttribute);

    template <typename K, typename V, typename C, typename A>
    void deserialize(const std::string& key, std::map<K, V*, C, A>& sMap,
                     const std::string& itemKey,
                     const std::string& comparisionAttribute = SerializeConstants::KeyAttribute);

    template <typename K, typename V, typename C, typename A>
    void deserialize(const std::string& key, std::map<K, std::unique_ptr<V>, C, A>& sMap,
                     const std::string& itemKey,
                     const std::string& comparisionAttribute = SerializeConstants::KeyAttribute);



    template <typename T, typename K>
    void deserialize(const std::string& key, ContainerWrapper<T, K>& container);


    // Specializations for chars
    void deserialize(const std::string& key, signed char& data, const bool asAttribute = false);
    void deserialize(const std::string& key, char& data, const bool asAttribute = false);
    void deserialize(const std::string& key, unsigned char& data, const bool asAttribute = false);

    // integers, reals, strings
    template <typename T, typename std::enable_if<std::is_integral<T>::value ||
                                                      std::is_floating_point<T>::value ||
                                                      util::is_string<T>::value,
                                                  int>::type = 0>
    void deserialize(const std::string& key, T& data, const bool asAttribute = false);

    // Enum types
    template <typename T, typename std::enable_if<std::is_enum<T>::value, int>::type = 0>
    void deserialize(const std::string& key, T& data, const bool asAttribute = false);

    // glm vector types
    template <typename Vec, typename std::enable_if<util::rank<Vec>::value == 1, int>::type = 0>
    void deserialize(const std::string& key, Vec& data);

    // glm matrix types
    template <typename Mat, typename std::enable_if<util::rank<Mat>::value == 2, int>::type = 0>
    void deserialize(const std::string& key, Mat& data);

    /**
     * \brief  Deserialize any Serializable object
     */
    void deserialize(const std::string& key, Serializable& sObj);
    /**
     * \brief  Deserialize pointer data of type T, which is of type
     *         serializeble object or primitive data
     */
    template <class T>
    void deserialize(const std::string& key, T*& data);

    void convertVersion(VersionConverter* converter);
    
    /**
     * \brief For allocating objects such as processors, properties.. using registered factories.
     *
     * @param className is used by registered factories to allocate the required object.
     * @return T* nullptr if allocation fails or className does not exist in any factories.
     */
    template <typename T>
    T* getRegisteredType(const std::string& className);

    /**
     * \brief For allocating objects that do not belong to any registered factories.
     *
     * @return T* Pointer to object of type T.
     */
    template <typename T>
    T* getNonRegisteredType();

    friend class NodeSwitch;
    template <typename T, typename K>
    friend class ContainerWrapper;

private:
    void registerFactories(InviwoApplication* app);

    void storeReferences(TxElement* node);

    void handleError(SerializationException&);

    std::vector<BaseDeserializationErrorHandler*> errorHandlers_;
    std::map<std::string, TxElement*> referenceLookup_;
    
    std::vector<FactoryBase*> registeredFactories_;
};

/**
 * \class ContainerWrapper
 */
template <typename T, typename K = std::string>
class ContainerWrapper {
public:
    using ItemAndCallback = std::pair<T&, std::function<void(T&)>>;
    using Getter = std::function<ItemAndCallback(K id, size_t ind)>;
    using IdentityGetter = std::function<K(TxElement* node)>;

    ContainerWrapper(std::string itemKey, Getter getItem) : itemKey_(itemKey), getItem_(getItem) {}

    virtual ~ContainerWrapper() = default;

    const std::string& getItemKey() const { return itemKey_; }

    void deserialize(Deserializer& d, TxElement* node, size_t ind) {
        auto itemAndCallback = getItem_(idGetter_(node), ind);
        try {
            d.deserialize(itemKey_, itemAndCallback.first);
            itemAndCallback.second(itemAndCallback.first);
        } catch (SerializationException& e) {
            d.handleError(e);
        }
    }

    void setIdentityGetter(IdentityGetter getter) { idGetter_ = getter; }

private:
    IdentityGetter idGetter_ = [](TxElement* node) {
        K key{};
        node->GetAttribute("identifier", &key, false);
        return key;
    };

    Getter getItem_;
    const std::string itemKey_;
};

namespace util {

template <typename T>
class IndexedDeserializer {
public:
    IndexedDeserializer(std::string key, std::string itemKey) : key_(key), itemKey_(itemKey) {}

    IndexedDeserializer<T>& setMakeNew(std::function<T()> makeNewItem) {
        makeNewItem_ = makeNewItem;
        return *this;
    }
    IndexedDeserializer<T>& onNew(std::function<void(T&)> onNewItem) {
        onNewItem_ = onNewItem;
        return *this;
    }
    IndexedDeserializer<T>& onRemove(std::function<bool(T&)> onRemoveItem) {
        onRemoveItem_ = onRemoveItem;
        return *this;
    }

    template <typename C>
    void operator()(Deserializer& d, C& container) {
        T tmp{};
        size_t count = 0;
        ContainerWrapper<T> cont(itemKey_, [&](std::string id, size_t ind) ->
                                 typename ContainerWrapper<T>::ItemAndCallback {
                                     ++count;
                                     if (ind < container.size()) {
                                         return {container[ind], [&](T& val) {}};
                                     } else {
                                         tmp = makeNewItem_();
                                         return {tmp, [&](T& val) { onNewItem_(val); }};
                                     }
                                 });

        d.deserialize(key_, cont);

        size_t n = 0;
        util::erase_remove_if(container, [&](T& item) {
            if (n < count) {
                ++n;
                return false;
            } else {
                ++n;
                return onRemoveItem_(item);
            }
        });
    }

private:
    std::function<T()> makeNewItem_ = []() -> T {
        throw Exception("MakeNewItem callback is not set!");
    };
    std::function<void(T&)> onNewItem_ = [](T&) {
        throw Exception("OnNewItem callback is not set!");
    };
    std::function<bool(T&)> onRemoveItem_ = [](T&) { return true; };

    std::string key_;
    std::string itemKey_;
};

template <typename K, typename T>
class IdentifiedDeserializer {
public:
    IdentifiedDeserializer(std::string key, std::string itemKey) : key_(key), itemKey_(itemKey) {}

    IdentifiedDeserializer<K, T>& setGetId(std::function<K(const T&)> getID) {
        getID_ = getID;
        return *this;
    }
    IdentifiedDeserializer<K, T>& setMakeNew(std::function<T()> makeNewItem) {
        makeNewItem_ = makeNewItem;
        return *this;
    }
    IdentifiedDeserializer<K, T>& onNew(std::function<void(T&)> onNewItem) {
        onNewItem_ = onNewItem;
        return *this;
    }
    IdentifiedDeserializer<K, T>& onRemove(std::function<void(const K&)> onRemoveItem) {
        onRemoveItem_ = onRemoveItem;
        return *this;
    }

    template <typename C>
    void operator()(Deserializer& d, C& container) {
        T tmp{};
        auto toRemove = util::transform(container, [&](const T& x)->K {return getID_(x);});
        ContainerWrapper<T, K> cont(
            itemKey_, [&](K id, size_t ind) -> typename ContainerWrapper<T, K>::ItemAndCallback {
                util::erase_remove(toRemove, id);
                auto it = util::find_if(container, [&](T& i) { return getID_(i) == id; });
                if (it != container.end()) {
                    return {*it, [&](T& val) {}};
                } else {
                    tmp = makeNewItem_();
                    return {tmp, [&](T& val) { onNewItem_(val); }};
                }
            });

        d.deserialize(key_, cont);
        for (auto& id : toRemove) onRemoveItem_(id);
    }

private:
    std::function<K(const T&)> getID_ = [](const T&) -> K { throw Exception("GetID callback is not set!"); };
    std::function<T()> makeNewItem_ = []() -> T {
        throw Exception("MakeNew callback is not set!");
    };
    std::function<void(T&)> onNewItem_ = [](T&) { throw Exception("OnNew callback is not set!"); };
    std::function<void(const K&)> onRemoveItem_ = [](const K&) {
        throw Exception("OnRemove callback is not set!");
    };

    std::string key_;
    std::string itemKey_;
};

template <typename K, typename T>
class MapDeserializer {
public:
    MapDeserializer(std::string key, std::string itemKey) : key_(key), itemKey_(itemKey) {}

    MapDeserializer<K, T>& setMakeNew(std::function<T()> makeNewItem) {
        makeNewItem_ = makeNewItem;
        return *this;
    }
    MapDeserializer<K, T>& onNew(std::function<void(const K&, T&)> onNewItem) {
        onNewItem_ = onNewItem;
        return *this;
    }
    MapDeserializer<K, T>& onRemove(std::function<void(const K&)> onRemoveItem) {
        onRemoveItem_ = onRemoveItem;
        return *this;
    }

    template <typename C>
    void operator()(Deserializer& d, C& container) {
        T tmp{};
        auto toRemove =
            util::transform(container, [](const std::pair<const K, T>& item) { return item.first; });
        ContainerWrapper<T, K> cont(
            itemKey_, [&](K id, size_t ind) -> typename ContainerWrapper<T, K>::ItemAndCallback {
                util::erase_remove(toRemove, id);
                auto it = container.find(id);
                if (it != container.end()) {
                    return {it->second, [&](T& val) {}};
                } else {
                    tmp = makeNewItem_();
                    return {tmp, [&, id](T& val) { onNewItem_(id, val); }};
                }
            });

        cont.setIdentityGetter([](TxElement* node) {
            K key{};
            node->GetAttribute(SerializeConstants::KeyAttribute, &key);
            return key;
        });

        d.deserialize(key_, cont);

        for (auto& id : toRemove) onRemoveItem_(id);
    }

private:
    std::function<T()> makeNewItem_ = []() -> T {
        throw Exception("MakeNew callback is not set!");
    };
    std::function<void(const K&, T&)> onNewItem_ = [](const K&, T&) {
        throw Exception("OnNew callback is not set!");
    };
    std::function<void(const K&)> onRemoveItem_ = [](const K&) {
        throw Exception("OnRemove callback is not set!");
    };

    std::string key_;
    std::string itemKey_;
};

}  // namespace


template <typename T>
T* Deserializer::getRegisteredType(const std::string& className) {
    for (auto base : registeredFactories_) {
        if (auto factory = dynamic_cast<Factory<T>*>(base)) {
            if (auto data = factory->create(className)) return data.release();
        }
    }
    return nullptr;
}

template <typename T>
T* Deserializer::getNonRegisteredType() {
    return util::defaultConstructType<T>();
}

template <typename T>
class DeserializationErrorHandle {
public:
    typedef void (T::*Callback)(SerializationException&);
    DeserializationErrorHandle(Deserializer&, std::string type, T* obj, Callback callback);
    virtual ~DeserializationErrorHandle();

private:
    Deserializer& d_;
};

template <typename T>
inviwo::DeserializationErrorHandle<T>::DeserializationErrorHandle(Deserializer& d,
                                                                  std::string type, T* obj,
                                                                  Callback callback)
    : d_(d) {
    d_.pushErrorHandler(new DeserializationErrorHandler<T>(type, obj, callback));
}

template <typename T>
inviwo::DeserializationErrorHandle<T>::~DeserializationErrorHandle() {
    delete d_.popErrorHandler();
}

// integers, reals, strings
template <typename T,
          typename std::enable_if<std::is_integral<T>::value || std::is_floating_point<T>::value ||
                                      util::is_string<T>::value,
                                  int>::type>
void Deserializer::deserialize(const std::string& key, T& data, const bool asAttribute) {
    if (asAttribute) {
        try {
            rootElement_->GetAttribute(key, &data);
        } catch (TxException&) {
        }
    } else {
        try {
            NodeSwitch ns(*this, key);
            rootElement_->GetAttribute(SerializeConstants::ContentAttribute, &data);
        } catch (TxException&) {
            try {
                NodeSwitch ns(*this, key, true);
                rootElement_->GetAttribute(SerializeConstants::ContentAttribute, &data);
            } catch (TxException&) {
            }
        }
    }
}

// enum types
template <typename T, typename std::enable_if<std::is_enum<T>::value, int>::type>
void Deserializer::deserialize(const std::string& key, T& data, const bool asAttribute) {
    using ET = typename std::underlying_type<T>::type;
    ET tmpdata{static_cast<ET>(data)};
    deserialize(key, tmpdata, asAttribute);
    data = static_cast<T>(tmpdata);
}

// glm vector types
template <typename Vec, typename std::enable_if<util::rank<Vec>::value == 1, int>::type>
void Deserializer::deserialize(const std::string& key, Vec& data) {
    try {
        NodeSwitch ns(*this, key);
        for (size_t i = 0; i < util::extent<Vec, 0>::value; ++i) {
            rootElement_->GetAttribute(SerializeConstants::VectorAttributes[i], &data[i]);
        }
    } catch (TxException&) {
    }
}

// glm matrix types
template <typename Mat, typename std::enable_if<util::rank<Mat>::value == 2, int>::type>
void Deserializer::deserialize(const std::string& key, Mat& data) {
    try {
        NodeSwitch ns(*this, key);
        for (size_t i = 0; i < util::extent<Mat, 0>::value; ++i) {
            deserialize("row" + toString(i), data[i]);
        }
    } catch (TxException&) {
    }
}

template <typename T>
void Deserializer::deserialize(const std::string& key, std::vector<T*>& vector,
                                  const std::string& itemKey) {
    try {
        NodeSwitch vectorNodeSwitch(*this, key);

        unsigned int i = 0;
        TxEIt child(itemKey);
        for (child = child.begin(rootElement_); child != child.end(); ++child) {
            // In the next deserialization call do net fetch the "child" since we are looping...
            // hence the "false" as the last arg.
            NodeSwitch elementNodeSwitch(*this, &(*child), false);

            if (vector.size() <= i) {
                T* item = nullptr;
                try {
                    deserialize(itemKey, item);
                    vector.push_back(item);
                } catch (SerializationException& e) {
                    delete item;
                    handleError(e);
                }
            } else {
                try {
                    deserialize(itemKey, vector[i]);
                } catch (SerializationException& e) {
                    handleError(e);
                }
            }
            i++;
        }
    } catch (TxException&) {
    }
}

template <typename T>
void Deserializer::deserialize(const std::string& key, std::vector<std::unique_ptr<T>>& vector,
                                  const std::string& itemKey) {
    try {
        NodeSwitch vectorNodeSwitch(*this, key);

        unsigned int i = 0;
        TxEIt child(itemKey);
        for (child = child.begin(rootElement_); child != child.end(); ++child) {
            // In the next deserialization call do net fetch the "child" since we are looping...
            // hence the "false" as the last arg.
            NodeSwitch elementNodeSwitch(*this, &(*child), false);

            if (vector.size() <= i) {
                T* item = nullptr;
                try {
                    deserialize(itemKey, item);
                    vector.emplace_back(item);
                } catch (SerializationException& e) {
                    delete item;
                    handleError(e);
                }
            } else {
                try {
                    auto ptr = vector[i].get();
                    deserialize(itemKey, ptr);
                } catch (SerializationException& e) {
                    handleError(e);
                }
            }
            i++;
        }
    } catch (TxException&) {
    }
}

template <typename T, typename C>
void Deserializer::deserialize(const std::string& key, std::vector<T*>& vector,
                                  const std::string& itemKey, C identifier) {
    try {
        NodeSwitch vectorNodeSwitch(*this, key);

        typename std::vector<T*>::iterator lastInsertion = vector.begin();

        TxEIt child(itemKey);
        for (child = child.begin(rootElement_); child != child.end(); ++child) {
            
            identifier.setKey(child.Get());
            auto it = std::find_if(vector.begin(), vector.end(), identifier);

            if (it != vector.end()) {  // There is a item in vector with same identifier as on disk
                NodeSwitch elementNodeSwitch(*this, &(*child), false);
                try {
                    deserialize(itemKey, *it);
                    lastInsertion = it;
                } catch (SerializationException& e) {
                    handleError(e);
                }
            } else {  // No item in vector matches item on disk, create a new one.
                T* item = nullptr;
                NodeSwitch elementNodeSwitch(*this, &(*child), false);
                try {
                    deserialize(itemKey, item);
                    // Insert new item after the previous item deserialized
                    lastInsertion = lastInsertion == vector.end() ? lastInsertion : ++lastInsertion;
                    lastInsertion = vector.insert(lastInsertion, item);
                } catch (SerializationException& e) {
                    delete item;
                    handleError(e);
                }
            }
        }
    } catch (TxException&) {
    }
}

template <typename T>
void Deserializer::deserialize(const std::string& key, std::vector<T>& vector,
                                  const std::string& itemKey) {
    try {
        NodeSwitch vectorNodeSwitch(*this, key);
        unsigned int i = 0;
        TxEIt child(itemKey);

        for (child = child.begin(rootElement_); child != child.end(); ++child) {
            // In the next deserialization call do net fetch the "child" since we are looping...
            // hence the "false" as the last arg.
            NodeSwitch elementNodeSwitch(*this, &(*child), false);

            try {
                if (vector.size() <= i) {
                    T item;
                    deserialize(itemKey, item);
                    vector.push_back(item);
                } else {
                    deserialize(itemKey, vector[i]);
                }
            } catch (SerializationException& e) {
                handleError(e);
            }
            i++;
        }
    } catch (TxException&) {
    }
}

template <typename T>
void Deserializer::deserialize(const std::string& key, std::list<T>& container,
                                  const std::string& itemKey) {
    try {
        NodeSwitch vectorNodeSwitch(*this, key);
        unsigned int i = 0;
        TxEIt child(itemKey);

        for (child = child.begin(rootElement_); child != child.end(); ++child) {
            // In the next deserialization call do net fetch the "child" since we are looping...
            // hence the "false" as the last arg.
            NodeSwitch elementNodeSwitch(*this, &(*child), false);

            try {
                if (container.size() <= i) {
                    T item;
                    deserialize(itemKey, item);
                    container.push_back(item);
                } else {
                    deserialize(itemKey, *std::next(container.begin(), i));
                }
            } catch (SerializationException& e) {
                handleError(e);
            }
            i++;
        }
    } catch (TxException&) {
    }
}

template <typename K, typename V, typename C, typename A>
void Deserializer::deserialize(const std::string& key, std::map<K, V, C, A>& map,
                                  const std::string& itemKey,
                                  const std::string& comparisionAttribute) {
    if (!isPrimitiveType(typeid(K)) || comparisionAttribute.empty())
        throw SerializationException("Error: map key has to be a primitive type", IvwContext);

    try {
        NodeSwitch mapNodeSwitch(*this, key);
        TxEIt child(itemKey);

        for (child = child.begin(rootElement_); child != child.end(); ++child) {
            // In the next deserialization call do net fetch the "child" since we are looping...
            // hence the "false" as the last arg.
            NodeSwitch elementNodeSwitch(*this, &(*child), false);
            K childkey;
            child->GetAttribute(comparisionAttribute, &childkey);

            V value;
            auto it = map.find(childkey);
            if(it != map.end()) value = it->second;
            try {
                deserialize(itemKey, value);
                map[childkey] = value;
            } catch (SerializationException& e) {
                handleError(e);
            }
        }
    } catch (TxException&) {
    }
}

// Pointer specialization
template <typename K, typename V, typename C, typename A>
void Deserializer::deserialize(const std::string& key, std::map<K, V*, C, A>& map,
                               const std::string& itemKey,
                               const std::string& comparisionAttribute) {
    if (!isPrimitiveType(typeid(K)) || comparisionAttribute.empty())
        throw SerializationException("Error: map key has to be a primitive type", IvwContext);

    try {
        NodeSwitch mapNodeSwitch(*this, key);
        TxEIt child(itemKey);

        for (child = child.begin(rootElement_); child != child.end(); ++child) {
            // In the next deserialization call do net fetch the "child" since we are looping...
            // hence the "false" as the last arg.
            NodeSwitch elementNodeSwitch(*this, &(*child), false);
            K childkey;
            child->GetAttribute(comparisionAttribute, &childkey);

            auto it = map.find(childkey);
            V*  value = (it != map.end() ? it->second : nullptr);

            try {
                deserialize(itemKey, value);
                map[childkey] = value;
            } catch (SerializationException& e) {
                handleError(e);
            }
        }
    } catch (TxException&) {
    }
}

// unique_ptr specialization
template <typename K, typename V, typename C, typename A>
void Deserializer::deserialize(const std::string& key, std::map<K, std::unique_ptr<V>, C, A>& map,
                               const std::string& itemKey,
                               const std::string& comparisionAttribute) {
    if (!isPrimitiveType(typeid(K)) || comparisionAttribute.empty())
        throw SerializationException("Error: map key has to be a primitive type", IvwContext);

    try {
        NodeSwitch mapNodeSwitch(*this, key);
        TxEIt child(itemKey);

        for (child = child.begin(rootElement_); child != child.end(); ++child) {
            // In the next deserialization call do net fetch the "child" since we are looping...
            // hence the "false" as the last arg.
            NodeSwitch elementNodeSwitch(*this, &(*child), false);
            K childkey;
            child->GetAttribute(comparisionAttribute, &childkey);

            auto it = map.find(childkey);
            if (it != map.end()){
                try {
                    deserialize(itemKey, it->second.get());
                 } catch (SerializationException& e) {
                    handleError(e);
                }
            } else {
                V* ptr = nullptr;
                try {      
                    deserialize(itemKey, ptr);
                    map.emplace(childkey, std::unique_ptr<V>(ptr));
                } catch (SerializationException& e) {
                    delete ptr;
                    handleError(e);
                }
            
            }
        }
    } catch (TxException&) {
    }
}


template <class T>
inline void Deserializer::deserialize(const std::string& key, T*& data) {
    TxElement* keyNode =
        retrieveChild_ ? rootElement_->FirstChildElement(key, false) : rootElement_;
    if (!keyNode) return;

    const std::string type_attr(keyNode->GetAttribute(SerializeConstants::TypeAttribute));
    const std::string ref_attr(keyNode->GetAttribute(SerializeConstants::RefAttribute));
    const std::string id_attr(keyNode->GetAttribute(SerializeConstants::IDAttribute));

    if (!data) {
        if (allowRef_ && !ref_attr.empty()) {
            // Has reference identifier, data should already be allocated and we only have to find
            // and set the pointer to it.
            data = static_cast<T*>(refDataContainer_.find(type_attr, ref_attr));
            if (!data) {
                std::map<std::string, TxElement*>::iterator it = referenceLookup_.find(ref_attr);
                if (it != referenceLookup_.end()) {
                    NodeDebugger error(keyNode);
                    throw SerializationException(
                        "Reference to " + error[0].key + " not instantiated: \"" +
                            error[0].identifier + "\" of class \"" + error[0].type + "\" at line " +
                            toString(error[0].line),
                        IvwContext, error[0].key, error[0].type, error[0].identifier, it->second);
                } else {
                    throw SerializationException(
                        "Could not find reference to " + key + ": " + type_attr, IvwContext, key,
                        type_attr);
                }
            }
            return;

        } else if (!type_attr.empty()) {
            data = getRegisteredType<T>(type_attr);
            if (!data) {
                NodeDebugger error(keyNode);
                throw SerializationException(
                    "Could not create " + key + ": \"" + error[0].identifier + "\" of class \"" +
                        error[0].type + "\" at line: " + toString(error[0].line),
                    IvwContext, key, type_attr, error[0].identifier, keyNode);
            }

        } else {
            data = getNonRegisteredType<T>();
            if (!data) {
                NodeDebugger error(keyNode);
                throw SerializationException(
                    "Could not create " + key + ": \"" + error[0].identifier + "\" of class \"" +
                        error[0].type + "\" at line: " + toString(error[0].line),
                    IvwContext, key, type_attr, error[0].identifier, keyNode);
            }
        }
    }

    if (data) {
        deserialize(key, *data);
        if (allowRef_ && !id_attr.empty()) refDataContainer_.insert(data, keyNode);
    }
}

template <typename T, typename K>
void Deserializer::deserialize(const std::string& key, ContainerWrapper<T, K>& container) {
    try {
        NodeSwitch vectorNodeSwitch(*this, key);
        unsigned int i = 0;
        TxEIt child(container.getItemKey());

        for (child = child.begin(rootElement_); child != child.end(); ++child) {
            // In the next deserialization call do net fetch the "child" since we are looping...
            // hence the "false" as the last arg.
            NodeSwitch elementNodeSwitch(*this, &(*child), false);
            container.deserialize(*this, &(*child), i);
            i++;
        }
    } catch (TxException&) {
    }
}

}  // namespace
#endif